import argparse
import os.path
from pathlib import Path
import pickle
from tqdm import tqdm

def dict_train_test_split(complex_list,train_ids,val_ids,test_ids,del_train_ids,del_test_ids):
    d_train=[]
    d_val=[]
    complex_ids=set()
    for c in complex_list:
        if c['name'] in complex_ids:
            continue
        complex_ids.add(c['name'])
        if c['name_long'] in train_ids:
            d_train.append(c)
        elif c['name_long'] in val_ids:
            d_val.append(c)
        elif c['name_long'] in test_ids:
            pass
        elif c['name_long'] in del_train_ids:
            pass
        elif c['name_long'] in del_test_ids:
            pass
        else:
            print(f'Complex id {c["name_long"]} is not presented in train or test IDs')
    return d_train,d_val

def main(args, preprocessed_path, train_ids, val_ids,test_ids,del_train_ids,del_test_ids,indexes_of_cdrs):
    import json, time, os, sys, glob
    import shutil
    import warnings
    import numpy as np
    import torch
    from torch import optim
    from torch.utils.data import DataLoader
    import queue
    import copy
    import torch.nn as nn
    import torch.nn.functional as F
    import random
    import os.path
    import subprocess
    from concurrent.futures import ProcessPoolExecutor    
    from utils import worker_init_fn, get_pdbs, loader_pdb, build_training_clusters, PDB_dataset, StructureDataset, StructureLoader
    from model_utils import featurize, loss_smoothed, loss_nll, get_std_opt, ProteinMPNN

    with open(preprocessed_path, "rb") as f:
         d = pickle.load(f)
    scaler = torch.cuda.amp.GradScaler()
     
    device = torch.device("cuda:0" if (torch.cuda.is_available()) else "cpu")

    base_folder = time.strftime(args.path_for_outputs, time.localtime())
    
    df={'epoch': [], 'step':[], 'time': [], 'train': [], 'valid':[], 'train_acc': [], 'valid_acc': [],
            (0,0):[],(0,1):[],(0,2):[],
            (1,0):[],(1,1):[],(1,2):[],}
    if base_folder[-1] != '/':
        base_folder += '/'
    if not os.path.exists(base_folder):
        os.makedirs(base_folder)
    subfolders = ['model_weights']
    for subfolder in subfolders:
        if not os.path.exists(base_folder + subfolder):
            os.makedirs(base_folder + subfolder)

    PATH = args.previous_checkpoint

    logfile = base_folder + 'log.txt'
    if not PATH:
        with open(logfile, 'w') as f:
            f.write('Epoch\tTrain\tValidation\n')

    data_path = args.path_for_training_data
    params = {
        "LIST"    : f"{data_path}/list.csv", 
        "VAL"     : f"{data_path}/valid_clusters.txt",
        "TEST"    : f"{data_path}/test_clusters.txt",
        "DIR"     : f"{data_path}",
        "DATCUT"  : "2030-Jan-01",
        "RESCUT"  : args.rescut, #resolution cutoff for PDBs
        "HOMO"    : 0.70 #min seq.id. to detect homo chains
    }


    LOAD_PARAM = {'batch_size': 1,
                  'shuffle': True,
                  'pin_memory':False,
                  'num_workers': 4}

   
    if args.debug:
        args.num_examples_per_epoch = 50
        args.max_protein_length = 1000
        args.batch_size = 1000

    model = ProteinMPNN(node_features=args.hidden_dim, 
                        edge_features=args.hidden_dim, 
                        hidden_dim=args.hidden_dim, 
                        num_encoder_layers=args.num_encoder_layers, 
                        num_decoder_layers=args.num_encoder_layers, 
                        k_neighbors=args.num_neighbors, 
                        dropout=args.dropout, 
                        augment_eps=args.backbone_noise)
    model.to(device)
    if PATH:
        checkpoint = torch.load(PATH)
        # total_step = checkpoint['step'] #write total_step from the checkpoint
        # epoch = checkpoint['epoch'] #write epoch from the checkpoint
        total_step = 0
        epoch = 0
        model.load_state_dict(checkpoint['model_state_dict'])
    else:
        total_step = 0
        epoch = 0

    optimizer = get_std_opt(model.parameters(), args.hidden_dim, total_step)


    # if PATH:
    #     optimizer.optimizer.load_state_dict(checkpoint['optimizer_state_dict'])

    d_train,d_test=dict_train_test_split(d,train_ids,val_ids,test_ids,del_train_ids,del_test_ids)
    dataset_train = StructureDataset(d_train, truncate=None, max_length=args.max_protein_length) 
    dataset_valid = StructureDataset(d_test, truncate=None, max_length=args.max_protein_length)
    
    loader_train = StructureLoader(dataset_train, batch_size=args.batch_size)
    loader_valid = StructureLoader(dataset_valid, batch_size=args.batch_size)
    
    reload_c = 0 
    max_val_acc=0
    for e in range(args.num_epochs):
        t0 = time.time()
        e = epoch + e
        model.train()
        train_sum, train_weights = 0., 0.
        train_acc = 0.
        for _, batch in enumerate(loader_train):
            start_batch = time.time()
            X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, device,mode='valid',cdr_indexes=indexes_of_cdrs[0])
            elapsed_featurize = time.time() - start_batch
            optimizer.zero_grad()
            mask_for_loss = mask*chain_M
            
            if args.mixed_precision:
                with torch.cuda.amp.autocast():
                    log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding_all)
                    _, loss_av_smoothed = loss_smoothed(S, log_probs, mask_for_loss)
        
                scaler.scale(loss_av_smoothed).backward()
                    
                if args.gradient_norm > 0.0:
                    total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), args.gradient_norm)

                scaler.step(optimizer)
                scaler.update()
            else:
                log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding_all)
                _, loss_av_smoothed = loss_smoothed(S, log_probs, mask_for_loss)
                loss_av_smoothed.backward()

                if args.gradient_norm > 0.0:
                    total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), args.gradient_norm)

                optimizer.step()
            
            loss, loss_av, true_false = loss_nll(S, log_probs, mask_for_loss)
        
            train_sum += torch.sum(loss * mask_for_loss).cpu().data.numpy()
            train_acc += torch.sum(true_false * mask_for_loss).cpu().data.numpy()
            train_weights += torch.sum(mask_for_loss).cpu().data.numpy()

            total_step += 1
        dict_val={}
        
        
        for ind,jnd in indexes_of_cdrs:
            model.eval()
            with torch.no_grad():
                validation_sum, validation_weights = 0., 0.
                validation_acc = 0.
                for _, batch in enumerate(loader_valid):
                    X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, device,mode='valid',cdr_indexes=[ind,jnd])
                    log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding_all)
                    mask_for_loss = mask*chain_M
                    loss, loss_av, true_false = loss_nll(S, log_probs, mask_for_loss)
                    
                    validation_sum += torch.sum(loss * mask_for_loss).cpu().data.numpy()
                    validation_acc += torch.sum(true_false * mask_for_loss).cpu().data.numpy()
                    validation_weights += torch.sum(mask_for_loss).cpu().data.numpy()
            
            # train_loss = train_sum / train_weights
            # train_accuracy = train_acc / train_weights
            # train_perplexity = np.exp(train_loss)
            validation_loss = validation_sum / validation_weights
            validation_accuracy = validation_acc / validation_weights
            validation_perplexity = np.exp(validation_loss)
            
            #train_perplexity_ = np.format_float_positional(np.float32(train_perplexity), unique=False, precision=3)     
            validation_perplexity_ = np.format_float_positional(np.float32(validation_perplexity), unique=False, precision=3)
            #train_accuracy_ = np.format_float_positional(np.float32(train_accuracy), unique=False, precision=3)
            validation_accuracy_ = np.format_float_positional(np.float32(validation_accuracy), unique=False, precision=3)
            df[(ind,jnd)].append(validation_accuracy_)
        with torch.no_grad():
            validation_sum, validation_weights = 0., 0.
            validation_acc = 0.
            for _, batch in enumerate(loader_valid):
                X, S, mask, lengths, chain_M, residue_idx, mask_self, chain_encoding_all = featurize(batch, device,mode='valid',cdr_indexes=[ind,jnd])
                log_probs = model(X, S, mask, chain_M, residue_idx, chain_encoding_all)
                mask_for_loss = mask*chain_M
                loss, loss_av, true_false = loss_nll(S, log_probs, mask_for_loss)
                
                validation_sum += torch.sum(loss * mask_for_loss).cpu().data.numpy()
                validation_acc += torch.sum(true_false * mask_for_loss).cpu().data.numpy()
                validation_weights += torch.sum(mask_for_loss).cpu().data.numpy()

        validation_loss = validation_sum / validation_weights
        validation_accuracy = validation_acc / validation_weights
        validation_perplexity = np.exp(validation_loss)

        train_loss = train_sum / train_weights
        train_accuracy = train_acc / train_weights
        train_perplexity = np.exp(train_loss)
        
        train_perplexity_ = np.format_float_positional(np.float32(train_perplexity), unique=False, precision=3)     
        validation_perplexity_ = np.format_float_positional(np.float32(validation_perplexity), unique=False, precision=3)
        train_accuracy_ = np.format_float_positional(np.float32(train_accuracy), unique=False, precision=3)
        validation_accuracy_ = np.format_float_positional(np.float32(validation_accuracy), unique=False, precision=3)
        if validation_accuracy>max_val_acc:
            checkpoint_filename_best = base_folder+'model_weights/epoch_best.pt'.format(e+1, total_step)
            torch.save({
                        'epoch': e+1,
                        'step': total_step,
                        'num_edges' : args.num_neighbors,
                        'noise_level': args.backbone_noise,
                        'model_state_dict': model.state_dict(),
                        'optimizer_state_dict': optimizer.optimizer.state_dict(),
                        }, checkpoint_filename_best)
            max_val_acc=validation_accuracy
        t1 = time.time()
        dt = np.format_float_positional(np.float32(t1-t0), unique=False, precision=1) 
        with open(logfile, 'a') as f:
            f.write(f'epoch: {e+1}, step: {total_step}, time: {dt}, train: {train_perplexity_}, valid: {validation_perplexity_}, train_acc: {train_accuracy_}, valid_acc: {validation_accuracy_}\n')
        print(f'epoch: {e+1}, step: {total_step}, time: {dt}, train: {train_perplexity_}, valid: {validation_perplexity_}, train_acc: {train_accuracy_}, valid_acc: {validation_accuracy_}')
        for ind,jnd in indexes_of_cdrs:
            print(f'({ind},{jnd}): {df[(ind,jnd)][-1]}',end=' ')

        print()
        df['epoch'].append(e+1)
        df['step'].append(total_step) 
        df['time'].append(dt)
        df['train'].append(train_perplexity_)
        df['valid'].append(validation_perplexity_)
        df['train_acc'].append(train_accuracy_)
        df['valid_acc'].append(validation_accuracy_)

        
        checkpoint_filename_last = base_folder+'model_weights/epoch_last.pt'.format(e+1, total_step)
        torch.save({
                    'epoch': e+1,
                    'step': total_step,
                    'num_edges' : args.num_neighbors,
                    'noise_level': args.backbone_noise,
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.optimizer.state_dict(),
                    }, checkpoint_filename_last)

        if (e+1) % args.save_model_every_n_epochs == 0:
            checkpoint_filename = base_folder+'model_weights/epoch{}_step{}.pt'.format(e+1, total_step)
            torch.save({
                    'epoch': e+1,
                    'step': total_step,
                    'num_edges' : args.num_neighbors,
                    'noise_level': args.backbone_noise, 
                    'model_state_dict': model.state_dict(),
                    'optimizer_state_dict': optimizer.optimizer.state_dict(),
                    }, checkpoint_filename)
    
    with open(base_folder+'val_accuracy_CDRs_no_ckpt.pkl', 'wb') as fp:
        pickle.dump(df,fp)
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    argparser.add_argument("--path_for_training_data", type=str, default="my_path/pdb_2021aug02", help="path for loading training data") 
    argparser.add_argument("--path_for_outputs", type=str, default="./exp_020", help="path for logs and model weights")
    argparser.add_argument("--previous_checkpoint", type=str, default="", help="path for previous model weights, e.g. file.pt")
    argparser.add_argument("--num_epochs", type=int, default=200, help="number of epochs to train for")
    argparser.add_argument("--save_model_every_n_epochs", type=int, default=10, help="save model weights every n epochs")
    argparser.add_argument("--reload_data_every_n_epochs", type=int, default=2, help="reload training data every n epochs")
    argparser.add_argument("--num_examples_per_epoch", type=int, default=1000000, help="number of training example to load for one epoch")
    argparser.add_argument("--batch_size", type=int, default=10000, help="number of tokens for one batch")
    argparser.add_argument("--max_protein_length", type=int, default=10000, help="maximum length of the protein complext")
    argparser.add_argument("--hidden_dim", type=int, default=128, help="hidden model dimension")
    argparser.add_argument("--num_encoder_layers", type=int, default=3, help="number of encoder layers") 
    argparser.add_argument("--num_decoder_layers", type=int, default=3, help="number of decoder layers")
    argparser.add_argument("--num_neighbors", type=int, default=48, help="number of neighbors for the sparse graph")   
    argparser.add_argument("--dropout", type=float, default=0.1, help="dropout level; 0.0 means no dropout")
    argparser.add_argument("--backbone_noise", type=float, default=0.2, help="amount of noise added to backbone during training")   
    argparser.add_argument("--rescut", type=float, default=3.5, help="PDB resolution cutoff")
    argparser.add_argument("--debug", type=bool, default=False, help="minimal data loading for debugging")
    argparser.add_argument("--gradient_norm", type=float, default=-1.0, help="clip gradient norm, set to negative to omit clipping")
    argparser.add_argument("--mixed_precision", type=bool, default=True, help="train with mixed precision")
    
    argparser.add_argument("--preprocessed_path", type=Path, default="/mnt/proteinmpnn/ProteinMPNN_preprocessed_chothia_proteinlib_logging.pickle")
    argparser.add_argument("--regions", type=str, default="H3")
    argparser.add_argument("--comment", type=str, default="")

    args = argparser.parse_args() 
    preprocessed_path = Path(args.preprocessed_path)
    regions=args.regions
    args.path_for_outputs=f'exp_{regions}{args.comment}'
    train_file=Path(f'train_val_test_{regions}/train_renamed_clusterRes_0.5_DB_CDR_{regions}.fasta_cluster.txt')
    val_file=Path(f'train_val_test_{regions}/val_renamed_clusterRes_0.5_DB_CDR_{regions}.fasta_cluster.txt')
    test_file=Path(f'train_val_test_{regions}/test_renamed_clusterRes_0.5_DB_CDR_{regions}.fasta_cluster.tsv')
    del_train_file=Path(f'train_val_test_{regions}/deleted_train_and_val_renamed_clusterRes_0.5_DB_CDR_{regions}.fasta_cluster.tsv')
    del_test_file=Path(f'train_val_test_{regions}/deleted_train_and_val_renamed_clusterRes_0.5_DB_CDR_{regions}.fasta_cluster.tsv')
    

    train_ids=train_file.read_text().splitlines()
    val_ids=val_file.read_text().splitlines()
    test_ids=test_file.read_text().splitlines()
    del_train_ids=del_train_file.read_text().splitlines()
    del_test_ids=del_test_file.read_text().splitlines()
    l=['H1','H2','H3']
    indexes_of_cdrs=[(0,l.index(regions))]
    print(indexes_of_cdrs)
    main(args, preprocessed_path,train_ids,val_ids,test_ids,del_train_ids,del_test_ids,indexes_of_cdrs)   



# with open('data.pickle', 'rb') as f:
#     data = pickle.load(f)
# with open('val_accuracy_CDRs_no_ckpt.pkl', 'rb') as fp:
#     d=pickle.load(fp)