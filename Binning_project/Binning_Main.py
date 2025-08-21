#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import networkx as nx
import numpy as np
from Bio.SeqIO.FastaIO import SimpleFastaParser
from networkx.algorithms.community import louvain_communities
from sklearn.metrics.pairwise import cosine_similarity

def main():
    parser = argparse.ArgumentParser(description="Contig binning using k-mer similarity, coverage difference and community detection")
    parser.add_argument("fasta_file", help="Path to the FASTA file (e.g., contigs.fasta)")
    parser.add_argument("kmer_file", help="Path to the 4-mer reverse complement feature CSV file (e.g., 4mer_ReverseNor.csv)")
    parser.add_argument("-o", "--output", default="Metagenomic_Bins.csv", help="Name of the output CSV file (default: Bins.csv)")
    
    args = parser.parse_args()
    
   # --- FASTA
    contigs = []
    lengths = []
    with open(args.fasta_file) as fasta_file:
        for title, sequence in SimpleFastaParser(fasta_file):
            contigs.append(title.split(None, 1)[0])  # First word is ID
            lengths.append(len(sequence))

    dfLength = pd.DataFrame({'Contigs': contigs,'Length': lengths})
    df = dfLength.copy()

    # --- KMER CSV
    df1 = pd.read_csv(args.kmer_file, header=None)



## Computing contig similarity

    def process_group(df, df1, min_len, max_len, top_n, score_threshold,rou):
        df_combined = pd.concat([df, df1], axis=1)
        df_filtered = df_combined[(df_combined.Length >= min_len) & (df_combined.Length <= max_len)]
        features = df_filtered.drop(['Contigs', 'Length', 0], axis=1).round(rou)
        node_ids = df_filtered['Contigs']
    
        cos_sim_matrix = cosine_similarity(features)
    
        top_frame = []
        for i, row in enumerate(cos_sim_matrix):
            ix_top_n = np.argsort(-row)[:top_n]
            for j in ix_top_n:
                if i != j:
                    top_frame.append((node_ids.iloc[i], node_ids.iloc[j], row[j]))
    
        result = pd.DataFrame(top_frame, columns=['title1', 'title2', 'Score'])
        return result[result['Score'] >= score_threshold].sort_values(by='Score', ascending=False)


    df7 = process_group(df, df1, 9000, 1000000, 4, 0.85,3)
    df8 = process_group(df, df1, 3500, 25000, 4, 0.85,3)
    df9 = process_group(df, df1, 1000, 9000, 4, 0.85,3)
    
### Merge

    frames = [df7,df8,df9]
    df2 = pd.concat(frames)
    df2['title11'] = df2['title1']
    df2['title22'] = df2['title2']
    df2[['Node_source','length_S']] = df2['title11'].str.split('_length_',expand=True)
    df2[['Node_target','length_T']] = df2['title22'].str.split('_length_',expand=True)
    df2['Node_source'] = df2['Node_source'].str.replace(r'\D', '', regex=True)
    df2['Node_target'] = df2['Node_target'].str.replace(r'\D', '', regex=True)
    df2[['length_S','cov_S']] = df2['length_S'].str.split('_cov_',expand=True)
    df2[['length_T','cov_T']] = df2['length_T'].str.split('_cov_',expand=True)
    df3 = df2.drop(['title11','title22','Node_source','Node_target'], axis=1)
    df3[['cov_S','cov_T','length_S','length_T']] = df3[['cov_S','cov_T','length_S','length_T']].astype(float)

###  Difference Coverage

    df3['COV'] = df3.apply(lambda x: x['cov_S'] if x['cov_S'] <= x['cov_T'] else x['cov_T'], axis=1)
    df3['COVMul'] = (df3['COV']*1.5)
    df4 = df3.loc[ (df3['cov_S'] <= df3['COVMul']) & (df3['cov_T'] <= df3['COVMul'])]
    df5 = df4.drop(['COV','COVMul'], axis=1)
    df6 = df5.loc[ (df5['length_S'] != df5['length_T'])]

    G2 = nx.from_pandas_edgelist(df6, source='title1', target='title2',edge_attr=["length_S", "length_T"])

# Remove contigs < threshold

    for component in list(nx.connected_components(G2)):
      if len(component) <= 4:
         G2.remove_nodes_from(component)

    connected_components = list(nx.connected_components(G2))
    count_5_node_subgraphs = 0

    for component in connected_components:
      if len(component) < 7:
        subgraph = G2.subgraph(component)
        if all(
            d.get("length_S", float("inf")) < 3500 and d.get("length_T", float("inf")) < 3500
            for u, v, d in subgraph.edges(data=True)
        ):
            G2.remove_nodes_from(component)
            count_5_node_subgraphs += 1
            
# Recompute connected components
    connected_components = list(nx.connected_components(G2))
    Countcomponent = len(connected_components)

# Count how many subgraphs consist of exactly 5 nodes
    count_5_node_subgraphs = sum(1 for component in connected_components if len(component) == 5)

    for component in connected_components:
      subgraph = G2.subgraph(component)
      num_edges = subgraph.number_of_edges()

      if Countcomponent <= 10:
          if num_edges < 6:
              G2.remove_nodes_from(component)
      else:
          if count_5_node_subgraphs > 9:
              if num_edges < 5:
                  G2.remove_nodes_from(component)


    df10 = nx.to_pandas_edgelist(G2)

# Structuring complex networks

    G_Community_detec = nx.from_pandas_edgelist(df10, source='source', target='target')

# Clustering # Run Louvain community detection

    Community_louvain = louvain_communities(G_Community_detec,  seed= 42)
    DF_Community_louvain = pd.DataFrame(Community_louvain).T
    DF_Community_louvain1 = DF_Community_louvain.loc[:, DF_Community_louvain.count(axis='rows') > 4]
    DF_melted_numeric = pd.melt(DF_Community_louvain1, var_name='Bin', value_name='Contig')
    DF_louvain = DF_melted_numeric.dropna()
    DF_louvain.to_csv(args.output, index=False)
    
if __name__ == "__main__":
    main()
