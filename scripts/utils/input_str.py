import pandas as pd
import numpy as np

def unique_patients(samples):
    df_samples = pd.DataFrame({ "samples" : samples, 
                            "patient" :["-".join(i.split("-")[:3]) for i in samples] })
    df_samples["plate"] = [i.split("-")[5] for i in df_samples["samples"]]
    df_samples = df_samples.sort_values(by = "plate", ascending= False) \
                .drop_duplicates(subset = "patient", keep = "first")
    return df_samples

def pair_sample(type):
    meta_data = pd.read_csv("../processed_data/meta/" + type + 
                            "_meta_filtered.csv")
    if type == "CRC":
        df_bn = pd.concat([pd.read_csv("../COAD/filtered_STR/COAD_bn.csv"),
                                pd.read_csv("../READ/filtered_STR/READ_bn.csv")]).reset_index()
        df_stn = pd.concat([pd.read_csv("../COAD/filtered_STR/COAD_stn.csv"),
                                pd.read_csv("../READ/filtered_STR/READ_stn.csv")]).reset_index()
        df_pt = pd.concat([pd.read_csv("../COAD/filtered_STR/COAD_pt.csv"),
                                pd.read_csv("../READ/filtered_STR/READ_pt.csv")]).reset_index()
    else: 
        df_bn = pd.read_csv("../" + type + "/filtered_STR/" + type + "_bn.csv")
        df_stn = pd.read_csv("../" + type + "/filtered_STR/" + type + "_stn.csv")
        df_pt = pd.read_csv("../" + type + "/filtered_STR/" + type + "_pt.csv")

    df_bn = pd.concat([df_bn, df_stn]).reset_index()
    df_bn = df_bn.loc[df_bn["sample"].isin(meta_data["name"]),]
    df_pt = df_pt.loc[df_pt["sample"].isin(meta_data["name"]),]

    # select unique samples 
    unique_bn = unique_patients(df_bn["sample"].unique())
    unique_pt = unique_patients(df_pt["sample"].unique())
    merge_sample = pd.merge(unique_bn, unique_pt, on = "patient")
    df_bn = df_bn.loc[df_bn["sample"].isin(merge_sample["samples_x"])]
    df_pt = df_pt.loc[df_pt["sample"].isin(merge_sample["samples_y"])]

    df_bn["patient"] = ["-".join(i.split("-")[:3]) for i in df_bn["sample"]]
    df_pt["patient"] = ["-".join(i.split("-")[:3]) for i in df_pt["sample"]]
    pt_ns = pd.merge(df_pt, df_bn, on = ["patient", "tmp_id"], suffixes=('_t', '_n'))

    # keep autosomal STRs only 
    pt_ns["chr"] = [i.split("_")[0] for i in pt_ns["tmp_id"]]
    pt_ns = pt_ns.loc[~pt_ns.chr.isin(["chrM", "chrX", "chrY"]),:]
    
    return pt_ns

def tumor_solid_sample(type):
    meta_data = pd.read_csv("../processed_data/meta/" + type + 
                            "_meta_filtered.csv")
    if type == "CRC":
        df_pt = pd.concat([pd.read_csv("../COAD/filtered_STR/COAD_pt.csv"),
                                pd.read_csv("../READ/filtered_STR/READ_pt.csv")]).reset_index()
        df_stn = pd.concat([pd.read_csv("../COAD/filtered_STR/COAD_stn.csv"),
                                pd.read_csv("../READ/filtered_STR/READ_stn.csv")]).reset_index()
    else: 
        df_pt = pd.read_csv("../" + type + "/filtered_STR/" + type + "_pt.csv")
        df_stn = pd.read_csv("../" + type + "/filtered_STR/" + type + "_stn.csv")

    df_pt = df_pt.loc[df_pt["sample"].isin(meta_data["name"]),]
    df_stn = df_stn.loc[df_stn["sample"].isin(meta_data["name"])]

    # select unique samples 
    unique_pt = unique_patients(df_pt["sample"].unique())
    unique_stn = unique_patients(df_stn["sample"].unique())
    
    df_pt = df_pt.loc[df_pt["sample"].isin(unique_pt["samples"])]
    df_stn = df_stn.loc[df_stn["sample"].isin(unique_stn["samples"])]
    df_pt["patient"] = ["-".join(i.split("-")[:3]) for i in df_pt["sample"]]
    df_stn["patient"] = ["-".join(i.split("-")[:3]) for i in df_stn["sample"]]
    
    # keep autosomal STRs only 
    df_pt["chr"] = [i.split("_")[0] for i in df_pt["tmp_id"]]
    df_pt = df_pt.loc[~df_pt.chr.isin(["chrM", "chrX", "chrY"]),:]
    df_stn["chr"] = [i.split("_")[0] for i in df_stn["tmp_id"]]
    df_stn = df_stn.loc[~df_stn.chr.isin(["chrM", "chrX", "chrY"]),:]
    
    return df_pt, df_stn

if __name__ == "__main__":
    main()