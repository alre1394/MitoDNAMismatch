#Import libraries

import pandas as pd
import io
import requests
from Bio.KEGG import REST


#read input_file
df_genes=pd.read_csv('genes.tsv',sep='\t')
df_genes.head()

#function to convert ncbi id to kegg
def conv_ncbi_to_kegg(ncbi_id):
    conv_url="https://rest.kegg.jp/conv/genes/ncbi-geneid:"
    res=requests.get(conv_url+ncbi_id).text
    return res



#function to get all pathways for a specified gene
def get_gene_pathways(gene_id):
    main_url="https://rest.kegg.jp/link/pathway/"
    res=requests.get(main_url+gene_id).text
    
    return res

#get information for all genes in kegg human database
genes =pd.read_table( io.StringIO(REST.kegg_list("hsa").read()),header=None)

genes.head()
#specific all pathways for kegg human database
result=REST.kegg_list('pathway','hsa').read()
result=pd.read_table(io.StringIO(result), header=None)
result.rename({0:'ID',1:'description'},axis=1,inplace=True)

allentries=[]

#Subseting: get the gene list and get their associated functional pathways in kegg
for gene_id,gene_name in df_genes.values:
    
    final_id=gene_id.strip().split(':')[-1]
    kegg_id=conv_ncbi_to_kegg(final_id).strip('\n').split('\t')[-1]
    pathways=get_gene_pathways(kegg_id)
    #gene_description=genes[genes[0]==kegg_id]
    for line in pathways.strip('\n').split('\n'):
        pathwayid=line.split('\tpath:')[-1]
        #print(pathwayid)
        description=result[result['ID'].str.contains(pathwayid)]['description'].values[0]
        subentry=[gene_id,gene_name,kegg_id,pathwayid,description]
        allentries.append(subentry)
        
    #print(res)
#convert result to dataframe
df=pd.DataFrame(allentries,columns=['GeneID','GeneName','KEGGID','KEGGPathwayID','PathwayDescription'])

df.to_csv("pathways.csv")
