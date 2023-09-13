from pathlib import Path
import pandas as pd
from argparse import ArgumentParser

if __name__=='__main__':
    parser=ArgumentParser('Join tsvs')
    parser.add_argument('--tsvs',type=Path,nargs='+')
    parser.add_argument('--joined_csv',type=Path,default=Path('joined_clusters.tsv'))
    args=parser.parse_args()
    tsvs=args.tsvs
    joined_csv=args.joined_csv
    df_list=[pd.read_csv(t,sep='\t',header=None,names=[t.stem,'id']) for t in tsvs]
    df_joined=df_list[0].copy()
    for df in df_list[1:]:
        df_joined=df_joined.merge(df,on='id',how='outer',suffixes=('_x', '_y'))
    cols=['id']
    for col in df_joined.columns:
       if col!='id':
           cols.append(col)
    df_joined=df_joined[cols]
    if joined_csv.exists():
        raise ValueError(f"file with path {joined_csv} already exists, change output_file variable or delete file {joined_csv}")
    else:
        df_joined.to_csv(joined_csv,sep='\t',index=False)
