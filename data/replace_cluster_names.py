import pandas as pd
from pathlib import Path
from argparse import ArgumentParser

def df_to_dict(df):
    d={}
    for row in df.iterrows():
        k,v=row[1].values
        if k not in d:
            d[k]=[v]
        else:
            d[k].append(v)
    return d

def get_sorted_names(d):
    d_sort={}
    for k,v in d.items():
        d_sort[k]=min(v)
    return d_sort

def rename_clusters(df,d_sorted):
    df_new=pd.DataFrame()
    for row in df.iterrows():
        k,_=row[1].values
        #df_new=pd.concat([df_new,pd.DataFrame({'0':d_sorted[k],'1':row[1].values[1]})])
        df_new=pd.concat([df_new,pd.DataFrame({0:d_sorted[k],1:row[1].values[1]},index=[0])],ignore_index=True)
    return df_new

if __name__=='__main__':
    parser=ArgumentParser('Rename clusters')
    parser.add_argument('--cluster_tsv',type=Path)
    args=parser.parse_args()
    cluster_csv=args.cluster_tsv
    renamed_csv=Path(f'{str(cluster_csv.parent)}/renamed_{cluster_csv.name}')
    df=pd.read_csv(cluster_csv,sep="\t")
    d=df_to_dict(df)
    d_sorted=get_sorted_names(d)
    df_clusters_renamed=rename_clusters(df,d_sorted)
    assert df_clusters_renamed.shape[0]==df.shape[0]
    if renamed_csv.exists():
        raise ValueError(f"file with path {renamed_csv} already exists, change output_file variable or delete file {renamed_csv}")
    else:
        df_clusters_renamed.to_csv(renamed_csv,sep='\t',header=False,index=False)
