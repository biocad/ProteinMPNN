import pandas as pd
from pathlib import Path

def filter_sample(clusters_csv,benchmark_clusters,t):
    cl_bench=benchmark_clusters.loc[:,'cluster'].value_counts()
    new_bench=benchmark_clusters['id'].tolist()
    bench_deleted=[]
    cluster_deleted=pd.DataFrame()
    for _,row in benchmark_clusters.iterrows():
        if row['size']/cl_bench[row['cluster']]>t:
            new_bench.remove(row['id'])
            bench_deleted.append(row['id'])
        else:
            cluster_deleted=pd.concat([cluster_deleted,clusters_csv[clusters_csv[clustering_regions].isin([row['cluster']])]])
            clusters_csv=clusters_csv[~clusters_csv[clustering_regions].isin([row['cluster']])]
            
    return clusters_csv,pd.Series(new_bench),cluster_deleted,pd.Series(bench_deleted)
     

if __name__=='__main__':
    benchmark_path=Path('rabd_benchmark.txt')
    clusters_path=Path('joined_clusters.tsv')
    clustering_regions='renamed_clusterRes_0.5_DB_CDR_H3.fasta_cluster'
    benchmark_csv=pd.read_csv(benchmark_path,header=None)
    clusters_csv=pd.read_csv(clusters_path,sep='\t')#.dropna()
    benchmark_clusters_csv=pd.DataFrame()
    l=[]
    s=0
    d=pd.DataFrame()
    for _,v in benchmark_csv.iterrows():
        for _,row in clusters_csv.iterrows():
                if v.values[0] in row['id']:
                    benchmark_clusters_csv.loc[row['id'],'id_complex']=v.values[0]
                    benchmark_clusters_csv.loc[row['id'],'id']=row['id']
                    benchmark_clusters_csv.loc[row['id'],'cluster']=row[clustering_regions]
                    vc=clusters_csv[clustering_regions].value_counts()
                    benchmark_clusters_csv.loc[row['id'],'size']=(row[clustering_regions]==clusters_csv[clustering_regions]).sum()
                    d.loc[row[clustering_regions],'s']=benchmark_clusters_csv.loc[row['id'],'size']
    clusters_csv[clustering_regions].isin(sorted(benchmark_clusters_csv.cluster.unique().tolist()))
    benchmark_clusters_csv['size']=benchmark_clusters_csv['size'].astype(int)
    benchmark_clusters_csv.sort_values('cluster').to_csv(f'benchmark_clusters_{clustering_regions}.tsv',sep='\t',index=False)

    

    clusters=clusters_csv[clustering_regions].unique()
    clusters_size=clusters_csv[clustering_regions].value_counts()
    in_benchmark=pd.DataFrame()

    for cl in clusters:
        ids_of_clusters=clusters_csv.loc[clusters_csv[clustering_regions]==cl]
        if set(ids_of_clusters[clustering_regions].tolist()).intersection(set(benchmark_clusters_csv['cluster'])):

            in_benchmark.loc[cl,'in_benchmark']=True#ids_of_clusters[clustering_regions].isin(benchmark_csv['cluster'])#.any()
            in_benchmark.loc[cl,'size']=benchmark_clusters_csv.loc[benchmark_clusters_csv['cluster']==cl,'size'].iloc[0]
        else:
            in_benchmark.loc[cl,'in_benchmark']=False
            in_benchmark.loc[cl,'size']=0
        #in_benchmark=pd.concat([in_benchmark,ser])
    in_benchmark['0']=in_benchmark.index
    in_benchmark['size']=in_benchmark['size'].astype(int)
    in_benchmark=in_benchmark[['0','in_benchmark','size']].dropna()
    in_benchmark.sort_values('0').to_csv(f'in_benchmark_{clustering_regions}.tsv',sep='\t',index=False)

    df=pd.DataFrame()
    l=[]
    # for _,row in in_bench.iterrows():
    #     #print(row['in_benchmark'])
    #     if row['in_benchmark']:
    #         l.append(row[0])
    ids_bench=in_benchmark.dropna().loc[in_benchmark['in_benchmark'],'0'].to_list()
    
    df_del=clusters_csv[~clusters_csv[clustering_regions].isin(ids_bench)]
    # for _,row in clusters_csv.iterrows():
    #     if row[clustering_regions] not in ids_bench:
    #         df=pd.concat([df,row.to_list()],ignore_index=True)

    t=30
    train,test,del_train,del_test=filter_sample(clusters_csv,benchmark_clusters_csv,t)
    train[['id',clustering_regions]].to_csv(f'train_and_val_{clustering_regions}.tsv',sep='\t',index=False)
    test.to_csv(f'test_{clustering_regions}.tsv',sep='\t',index=False,header=False)
    del_train[['id',clustering_regions]].to_csv(f'deleted_train_and_val_{clustering_regions}.tsv',sep='\t',index=False)
    del_test.to_csv(f'deleted_test_{clustering_regions}.tsv',sep='\t',index=False,header=False)
    #dropped_clusters.to_csv('dropped.tsv',index=False)
    #df_del.to_csv(f'filtered_sample_{clustering_regions}.tsv',sep='\t',index=False)
    # test
    s=set(train['id'].values.tolist())|set(test.values.tolist())|set(del_train['id'].values.tolist())|set(del_test.values.tolist())
    s=set(s)
    print(len(s))
    assert len(s)==pd.read_csv(clusters_path,sep='\t').shape[0]