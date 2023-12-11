import pandas as pd
from pathlib import Path

if __name__=='__main__':
    benchmark_path=Path('benchmark_clusters.tsv')
    clusters_path=Path('/mnt/sabdab/clusters.tsv')
    
    clusters_csv=pd.read_csv(clusters_path,sep='\t')#.dropna()
    benchmark_csv=pd.read_csv(benchmark_path,sep='\t')

    clusters=clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].unique()
    clusters_size=clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].value_counts()
    in_benchmark=pd.DataFrame()

    for cl in clusters:
        ids_of_clusters=clusters_csv.loc[clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster']==cl]
        if set(ids_of_clusters['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].tolist()).intersection(set(benchmark_csv['cluster'])):

            in_benchmark.loc[cl,'in_benchmark']=True#ids_of_clusters['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].isin(benchmark_csv['cluster'])#.any()
            in_benchmark.loc[cl,'size']=benchmark_csv.loc[benchmark_csv['cluster']==cl,'size'].iloc[0]
        else:
            in_benchmark.loc[cl,'in_benchmark']=False
            in_benchmark.loc[cl,'size']=0
        #in_benchmark=pd.concat([in_benchmark,ser])
    in_benchmark['0']=in_benchmark.index
    in_benchmark['size']=in_benchmark['size'].astype(int)
    in_benchmark=in_benchmark[['0','in_benchmark','size']]
    in_benchmark.sort_values('0').to_csv('in_benchmark.tsv',sep='\t',index=False)