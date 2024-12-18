import pandas as pd

## process STRs
def process_STR(str_path, stats_path, str_panel):
    
    name_to_sample = pd.read_csv(stats_path)[["sample", "file_name"]]
    df_str = pd.read_csv(str_path)
    df_str = df_str.merge(name_to_sample, left_on = "patient", right_on = "file_name", how = "left")
    df_str.drop(columns = ["patient", "file_name"], inplace = True)
    df_str.columns = [i.strip() for i in df_str.columns]
    
    df_str["tmp_id"] = df_str["chr"].str.cat(df_str["start"].astype("str"), sep = "_")
    
    df_to_dupe = (
        df_str[["sample", "tmp_id", "alt"]]
            .groupby(["sample", "tmp_id"])
            .filter(lambda x: len(x) == 1))
    
    df_wide = (
            pd.concat([df_to_dupe, df_str[["sample", "tmp_id", "alt"]]])
                .sort_values(by=["sample", "tmp_id"])
                .reset_index(drop=True))
    df_wide["allele"] = [f'allele_a', f'allele_b'] * int((df_wide.shape[0] / 2))
    
    df_wide = (df_wide
            .pivot(index=["sample", "tmp_id"], columns="allele", values="alt")
            .reset_index()
            .merge(df_str.drop(["sample", "alt"], axis=1).drop_duplicates(), how="left", on="tmp_id"))
    
    df_wide = df_wide.merge(str_panel, how = "inner", on = "tmp_id") \
                    .query("in_segdup == False and neighbour_type == 'no_neighbour'")
    df_out = df_wide[["sample", "tmp_id", "allele_a", "allele_b", "end_x", "period_x", "ref_x"]]
    df_out.columns = ["sample", "tmp_id", "allele_a", "allele_b", "end", "period", "ref"]
    
    return df_out

def main():
    cancer_type = "STAD"
    sample_type = "Blood_Derived_Normal"
    abbr = "bn"
    str_tral = pd.read_csv(".././STR_panel/tral_and_perf_panel_meta_info_updated.tsv", sep = "\t", low_memory = False)
    filtered_df = process_STR("./" + cancer_type + "/STR_" + cancer_type + "_" + sample_type + ".csv", 
                              "./" + cancer_type + "/stats_" + cancer_type + "_" + sample_type + ".csv", 
                              str_tral)
    filtered_df.to_csv("./" + cancer_type + "/filtered_STR/" + cancer_type + "_" + abbr + ".csv", index = False)

if __name__ == "__main__":
    main()
    
    