import pandas as pd
from pathlib import Path
from tqdm import tqdm

def filter_sample(clusters_csv,benchmark_clusters,t):
    #ids_bench=in_benchmark.dropna().loc[in_benchmark['in_benchmark'],'0'].to_list()
    cl_bench=benchmark_clusters.loc[:,'cluster'].value_counts()
    new_bench=benchmark_clusters['id'].tolist()
    for _,row in benchmark_clusters.iterrows():
        if row['size']/cl_bench[row['cluster']]>t:
            new_bench.remove(row['id'])
        else:
            clusters_csv=clusters_csv[~clusters_csv[clustering_regions].isin([row['cluster']])]
    return clusters_csv,pd.Series(new_bench)
     

if __name__=='__main__':
    # %%
    import pandas as pd
    from pathlib import Path
    from tqdm import tqdm

    def filter_sample(clusters_csv,benchmark_clusters,t):
        #ids_bench=in_benchmark.dropna().loc[in_benchmark['in_benchmark'],'0'].to_list()
        cl_bench=benchmark_clusters.loc[:,'cluster'].value_counts()
        new_bench=benchmark_clusters['id'].tolist()
        for _,row in benchmark_clusters.iterrows():
            if row['size']/cl_bench[row['cluster']]>t:
                new_bench.remove(row['id'])
            else:
                clusters_csv=clusters_csv[~clusters_csv[clustering_regions].isin([row['cluster']])]
        return clusters_csv,pd.Series(new_bench)

    clusters_path=Path('/mnt/sabdab/clusters_with_1w72.tsv')
    clusters_csv=pd.read_csv(clusters_path,sep='\t')#.dropna()

    clustering_regions='renamed_clusterRes_0.5_DB_CDR_H1_CDR_H2_CDR_H3_CDR_L1_CDR_L2_CDR_L3.fasta_cluster'
    clustering_regions='renamed_clusterRes_0.5_DB_CDR_H3.fasta_cluster'
    benchmark_path=f'/data/user/shapoval/ProteinMPNN/benchmark_clusters_{clustering_regions}.tsv'
    benchmark_clusters_csv=pd.read_csv(benchmark_path,sep='\t')
    
    T=range(200)
    train_size,test_size=[],[]
    for t in tqdm(T):
        train,test=filter_sample(clusters_csv,benchmark_clusters_csv,t)
        train_size.append(train.shape[0])
        test_size.append(test.shape[0])
    
    import matplotlib.pyplot as plt

    plt.plot(train_size,test_size)
    #plt.title(clustering_regions[19:67])
    plt.title(clustering_regions[19:])
    plt.grid()
    plt.xlabel('train size')
    plt.ylabel('test size')
    plt.show()
