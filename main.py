import scanpy as sc
import numpy as np
import json, os, subprocess
import argparse

parser = argparse.ArgumentParser(description="Process")
parser.add_argument("--dataset_path", type=str, help="Path of scRNAseq data",default="./annotated_seurat_filtered.rds")
parser.add_argument("--panel_num", type=int, help="Number of genes in panel", default=200)
parser.add_argument("--search_space_path", type=str, help="Path of searching gene set, needed to be a json file", default='./search.json')
parser.add_argument("--objmode", choices=["overall","cts","corr","pathway","spatial","tv"], required=True, help="Mode of chosing objective function")

args = parser.parse_args()

objmode = args.objmode

panel_num = args.panel_num
dataset_path = args.dataset_path
search_space_path= args.search_space_path
print(objmode)
panel_contains_path='panel_contains_'+objmode
panel_score_path='panel_score_'+objmode
os.makedirs(panel_contains_path, exist_ok=True)
os.makedirs(panel_score_path, exist_ok=True)

with open(search_space_path,'r') as f:
    all_genes=json.load(f)  
    print(len(all_genes))  
search_space = {}
for i in range(panel_num):
    search_space[str(i+1)]=all_genes

trial_command = 'Rscript'
script_path = 'scores.R'
max_trial_number = 2000
trial_concurrency = 10

from evo import *
evo_a=Evolution(panel_score_path, optimize_mode='maximize', population_size=50, objmode=objmode)
evo_a.update_search_space(search_space)

for i in range(max_trial_number // trial_concurrency):
    para_id_list=list(range(i*trial_concurrency+1,(i+1)*trial_concurrency+1))
    print(para_id_list)
    current_panels=evo_a.generate_multiple_panels(para_id_list)
    processes=[]
    for j in para_id_list:
        with open(os.path.join(panel_contains_path, str(j)+'.json'),'w') as f:
            json.dump(current_panels[j%trial_concurrency-1], f)
        
        panel_txt_path=os.path.join(panel_contains_path, str(j)+'.txt')
        with open(panel_txt_path, "w") as file:
            for gene_name in list(current_panels[j%trial_concurrency-1].values()):
                file.write(gene_name + "\n")
                
        cur_final_score_path=os.path.join(panel_score_path, str(j)+'_scores.csv')
        process = subprocess.Popen([trial_command, script_path, dataset_path, panel_txt_path, cur_final_score_path])
        processes.append(process)
    for process in processes:
        process.wait()
    for j in para_id_list:
        sucess_ornot, reward =evo_a.receive_trial_result(j)
        with open(os.path.join(panel_contains_path, str(j)),'w') as f:
            f.write(str(reward))
        evo_a.trial_end(j, sucess_ornot)