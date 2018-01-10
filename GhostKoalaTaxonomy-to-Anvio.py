#!/usr/bin/env python
import pandas as pd
import sys
GK_taxonomy = pd.read_table(sys.argv[1],header=None).replace({'genecall_': ''}, regex=True)
GK_taxonomy[0] = GK_taxonomy[0].map(lambda x: x.lstrip('user:'))
GK_taxonomy.columns=['gene_callers_id','acceession','t_phylum','t_class','t_genus','KeggGeneId','GHOSTX_Score']
extracted = GK_taxonomy.filter(['gene_callers_id','t_phylum','t_class','t_genus']).reset_index().set_index('gene_callers_id')
del extracted['index']
empty_row = [' ']*len(extracted['t_phylum'].tolist())
extracted.insert(2,'t_order',empty_row)
extracted.insert(3,'t_family',empty_row)
extracted.insert(5,'t_species',empty_row)
extracted.to_csv(sys.argv[2],sep='\t')
