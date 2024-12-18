import pandas as pd
import numpy as np
import glob

def combine_gene_file(folder):
    #folder = "/Users/bgadmin/Documents/gene_exp/CRC_gene/**/*.tsv"
    files = glob.glob(folder)
    first_file = files[0]
    single_sample = pd.read_csv(first_file, sep = "\t", comment = "#")
    first_gene = single_sample.iloc[4:,:].query(" gene_type == 'protein_coding' ")[["gene_id", "tpm_unstranded"]]
    file_name = first_file.split("/")[-1]
    first_gene.columns = ["gene_id", file_name]

    for file in files[1:]:
        single_sample = pd.read_csv(file, sep = "\t", comment = "#")
        single_gene = single_sample.iloc[4:,:].query(" gene_type == 'protein_coding' ")[["gene_id", "tpm_unstranded"]]
        file_name = file.split("/")[-1]
        single_gene.columns = ["gene_id", file_name]
        first_gene = pd.merge(first_gene, single_gene, on = "gene_id")
    first_gene.to_csv("CRC_gene_tpm.tsv", sep = "\t", index = False)

def file_to_id(file, file2id):
    gene_file = pd.read_csv(file, sep = "\t")
    file_names = pd.DataFrame({"file_name" : gene_file.columns[1:]})
    filetoid = pd.read_csv(file2id, sep = "\t")
    new_column = file_names.merge(filetoid, on = "file_name",how = "left")["sample_id"].values
    gene_file.columns = np.concatenate((["gene_id"], new_column))
    gene_file.to_csv(file, index = False, sep = "\t")

def main():
    
    meta = pd.read_json("./gene_exp/metadata_stad.cart.2024-04-10.json")
    name_id = pd.DataFrame(columns = ["file_name", "sample_id"])
    for i in range(meta.shape[0]):
        sample = meta.iloc[i,:]
        file_name = sample["file_name"]
        sample_id = sample["associated_entities"][0]['entity_submitter_id']
        name_id = pd.concat([name_id,
                            pd.DataFrame({"file_name" : [file_name],
                                        "sample_id" : [sample_id]})])
    name_id.to_csv("./gene_exp/stad_file2id.txt", index = False, sep = "\t")
        
    
if __name__ == "__main__":
    main()