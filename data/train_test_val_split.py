
from argparse import ArgumentParser
from pathlib import Path
import random
import pandas as pd


def train_test_split(df,test_part=0.2,seed=42):
    random.seed(seed)
    train,test=[],[]
    df.columns=[0,1]
    clusters=df[1].value_counts()
    i=0
    while len(test)<test_part*df.shape[0] and i<len(clusters):
        p=random.random()
        I=clusters.index[i]==df[1]
        if p>test_part:
            train+=df.loc[I,0].values.tolist()
        else:
            test+=df.loc[I,0].values.tolist()
        i+=1
    if i<len(clusters)-1:
        for j in range(i,len(clusters)):
            I=clusters.index[j]==df[1]
            train+=df.loc[I,0].values.tolist()

    return train, test


if __name__=='__main__':
    parser=ArgumentParser('Rename clusters')
    parser.add_argument('--cluster',default='clusterRes_0.5_DB_CDR_H3.fasta_cluster',type=str)
    parser.add_argument('--sample_dir',default=Path('train_val_test'))
    args=parser.parse_args()
    
    cluster=args.cluster
    sample_dir=Path(args.sample_dir)

    summary_csv=sample_dir/f'train_and_val_{cluster}.tsv'
    df = pd.read_csv(summary_csv,sep="\t").dropna()
    train,test=train_test_split(df)
    print(len(train)/df.shape[0])
    assert len(train)+len(test)==df.shape[0]
    pd.Series(train).to_csv(sample_dir/f'train_{cluster}.txt',index=False,header=False)
    pd.Series(test).to_csv(sample_dir/f'val_{cluster}.txt',index=False,header=False)