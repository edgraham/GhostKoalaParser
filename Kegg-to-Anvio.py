import pandas as pd
import argparse
"""
Parses annotation results from KEGG and optionally will pull in results from interproscan
Assumed interproscan was ran using the following flags: -f tsv --goterms --iprlookup --pathways
"""
parser = argparse.ArgumentParser(description='Combines annotation Data for input to anvio')
parser.add_argument('--KeggDB', 
                    help='identify the Kegg Orthology file (modified from htext using given bash script)')
parser.add_argument('-i', 
                    help='specify the file containing GhostKoala Results')
parser.add_argument('--interproscan',
                    help='interproscan results')
parser.add_argument('-o',
                    help='Specify an output file')                    
args = parser.parse_args()
arg_dict=vars(args)
keggortho_database = arg_dict['KeggDB']
output = arg_dict['o']
GK_results = arg_dict['i']
##Read In KO_Orthology file and format for downstream analysis##
x= pd.read_table(keggortho_database,header=None) 
y =pd.DataFrame(x[3].str.split(' ',1).tolist(),columns=['accession','description'])
xy = pd.concat([x,y],axis=1).drop(3,1).set_index('accession')
xy.columns= ["Category1","Category2","Category3","description"]
xy.to_csv("KeggOrthology_Table1.txt",encoding='utf-8')
keggAnnotation = pd.read_table(GK_results,header=None)
keggAnnotation.columns= ["gene_caller_id","accession"]
keggAnnotation=keggAnnotation.dropna().set_index("accession")
merged = keggAnnotation.join(xy)
merged_reduced = merged.drop_duplicates(subset='gene_caller_id', keep="last")
extracted = merged_reduced.filter(['gene_callers_id','function','accession']).reset_index().set_index('gene_caller_id')
e_value = [0]*len(extracted['accession'].tolist())
source = ['KeggGhostKoala']*len(extracted['accession'].tolist())
extracted.insert(0,'source',source)
extracted.insert(3,'e_value',e_value)
extracted = extracted.reset_index()
if arg_dict["interproscan"] is not None:
    interpro = pd.read_table(arg_dict["interproscan"],header=None)
    interpro.columns=["gene_caller_id","MD5","Length","source","accession","description","start_loc","stop_loc","e_value","status","date","InterProAccession","InterProDescription","GOAnnotations","Pathway"]
    InterProExtracted = interpro.filter(["gene_caller_id","source","accession","description","e_value"])
    KEGG_InterPro_Combined = pd.concat([extracted,InterProExtracted]).set_index('gene_caller_id')
    KEGG_InterPro_Combined.to_csv(output,sep='\t')
else:
    extracted.to_csv(output,sep='\t')