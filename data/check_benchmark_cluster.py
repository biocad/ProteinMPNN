import pandas as pd
from pathlib import Path

if __name__=='__main__':
    benchmark_path=Path('rabd_benchmark.txt')
    clusters_path=Path('/mnt/sabdab/clusters.tsv')
    #clustering_regions='renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'

    benchmark_csv=pd.read_csv(benchmark_path,header=None)
    clusters_csv=pd.read_csv(clusters_path,sep='\t')
    benchmark_clusters_csv=pd.DataFrame()
    l=[]
    s=0
    d=pd.DataFrame()
    for _,v in benchmark_csv.iterrows():
        for _,row in clusters_csv.iterrows():
                if v.values[0] in row['id']:
                    benchmark_clusters_csv.loc[row['id'],'id_complex']=v.values[0]
                    benchmark_clusters_csv.loc[row['id'],'id']=row['id']
                    benchmark_clusters_csv.loc[row['id'],'cluster']=row['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster']
                    vc=clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].value_counts()
                    benchmark_clusters_csv.loc[row['id'],'size']=(row['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster']==clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster']).sum()
                    d.loc[row['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'],'s']=benchmark_clusters_csv.loc[row['id'],'size']
    clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].isin(sorted(benchmark_clusters_csv.cluster.unique().tolist()))
    benchmark_clusters_csv['size']=benchmark_clusters_csv['size'].astype(int)
    benchmark_clusters_csv.sort_values('cluster').to_csv('benchmark_clusters.tsv',sep='\t',index=False)


