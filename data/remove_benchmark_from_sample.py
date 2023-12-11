import pandas as pd

if __name__=='__main__':
    in_bench=pd.read_csv('in_benchmark.tsv',sep='\t').dropna()
    clusters_csv=pd.read_csv('/mnt/sabdab/clusters.tsv',sep='\t').dropna()
    df=pd.DataFrame()
    l=[]
    # for _,row in in_bench.iterrows():
    #     #print(row['in_benchmark'])
    #     if row['in_benchmark']:
    #         l.append(row[0])
    ids_bench=in_bench.loc[in_bench['in_benchmark'],'0'].to_list()
    
    df_del=clusters_csv[~clusters_csv['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'].isin(ids_bench)]
    # for _,row in clusters_csv.iterrows():
    #     if row['renamed_clusterRes_0.5_DB_CDR_H3_CDR_H2_CDR_H1_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'] not in ids_bench:
    #         df=pd.concat([df,row.to_list()],ignore_index=True)
    df_del.to_csv('filtered_sample.tsv',sep='\t',index=False)