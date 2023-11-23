from pathlib import Path
import pandas as pd
from argparse import ArgumentParser

if __name__=='__main__':
    parser=ArgumentParser('Join tsvs')
    parser.add_argument('--tsvs',type=Path,nargs='+')
    parser.add_argument('--joined_tsv',type=Path,default=Path('MMseq/joined_clusters.tsv'))
    args=parser.parse_args()
    tsvs=args.tsvs
    joined_tsv=args.joined_tsv
    df_list=[pd.read_csv(t,sep='\t',header=None,names=[t.stem,'id']) for t in tsvs]
    df_joined=df_list[0].copy()
    for df in df_list[1:]:
        df_joined=df_joined.merge(df,on='id',how='outer',suffixes=('_x', '_y'))
    cols=['id']
    for col in df_joined.columns:
       if col!='id':
           cols.append(col)
    df_joined=df_joined[cols]
    if joined_tsv.exists():
        raise ValueError(f"file with path {joined_tsv} already exists, change output_file variable or delete file {joined_csv}")
    else:
        df_joined.to_csv(joined_tsv,sep='\t',index=False)
