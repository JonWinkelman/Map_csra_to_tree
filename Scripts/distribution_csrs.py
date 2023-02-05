#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 15 09:56:51 2022

@author: jonwinkelman
"""
import ete3
from ete3 import Tree
import math
from ete3 import NCBITaxa
import pandas as pd
import pandas
import random
import numpy as np
ncbi = NCBITaxa()
from Bio import Phyla

a=3
def convert_taxids_to_diff_rank():
    pass

def map_taxa_to_bac120(taxids):
    pass

df = pd.read_csv('/Volumes/JDW_ExtDrive/Lauren_Palmer/AstA_comparison_ornithineVsArginine/16s_geneTree/bac120/bac120_taxonomy.tsv', sep = '\t', header=None)
df.columns = ['seq_id', 'lineage']
df = df.set_index('seq_id')
t = Tree('/Volumes/JDW_ExtDrive/Lauren_Palmer/AstA_comparison_ornithineVsArginine/16s_geneTree/bac120/bac120.tree', format=1, quoted_node_names=True)
leaf_names = t.get_leaf_names()
df = df.loc[leaf_names,:]
a = [row.split(';') for row in df['lineage']]
df2 = pd.DataFrame(a)
df2.columns = ['domain', 'phylum','class','order','family','genus', 'species']
df2['seq_id'] = df.index
genus = [name[3:] for name in df2['genus']]
df2['genus'] = genus
genus_taxid_dict = ncbi.get_name_translator(genus)
genus_taxid = [genus_taxid_dict.get(name) for name in df2['genus']]
df2['genus_taxid'] =  [genus[0] if genus else None for genus in genus_taxid ]
 
csrs_df = pd.read_csv('csrs.tsv', sep='\t')
csr_genus_taxids = []
for taxon in csrs_df['Taxonomic lineage IDs']:
    if not math.isnan(taxon):
        lineage = ncbi.get_rank(ncbi.get_lineage(taxon))
        csr_genus_taxids.append([taxid for taxid, rank in lineage.items() if rank == 'genus'])
    else:
        csr_genus_taxids.append([])
csrs_df['genus_taxid'] = csr_genus_taxids
#remove list
for i, taxon in enumerate(csrs_df['genus_taxid']):
    if len(taxon) == 1:
        taxon = taxon[0]
    elif len(taxon) == 0:
        taxon = None
    else:
        print(taxon)
        taxon = None
    csrs_df.iloc[i,9] = taxon
not_nulls_filt = pd.notna(df2['genus_taxid'])
df2_filtered = df2.loc[not_nulls_filt,:]
unique_genuses = df2_filtered['genus_taxid'].unique()

acc_list = []
df2_i = df2.set_index('genus_taxid')
unique_genuses = [int(taxon) for taxon in unique_genuses]
for taxon in unique_genuses:
    a = df2_i.loc[taxon, 'seq_id']
    if type(a) == pandas.core.series.Series:
        a = list(a)
        b =[ele for ele in a if ele.startswith('RS_')]
        if len(b)<1:
            b = a
        samp = random.sample(b,1)[0]
        acc_list.append(samp)
    elif type(a) == str:
        acc_list.append(a)
        
df3 = df2.set_index('seq_id')
df3 = df3.loc[acc_list,:]      

bac120_genus_taxids = list(df3['genus_taxid'])
bac120_genus_taxids = [int(ele) for ele in bac120_genus_taxids]
b120_g_dict = {g_taxid:0 for g_taxid in bac120_genus_taxids}
for genus_taxid in csrs_df['genus_taxid']:
    if genus_taxid in b120_g_dict.keys():
        b120_g_dict[genus_taxid] = b120_g_dict[genus_taxid] + 1

t.prune(acc_list)
t.write(format=1, outfile = 'genus_bac120tree.nw')


df3 = df3.reset_index()
df3 = df3.set_index('genus_taxid')   
df4 = pd.DataFrame.from_dict(b120_g_dict, orient='index', columns=['counts'])
df5 = pd.merge(df3,df4, left_index=True, right_index=True)
df5.set_index('seq_id')['counts'].to_csv('genus-accession_counts.csv')