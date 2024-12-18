import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
import statsmodels.formula.api as sm
from statsmodels.stats import multitest

def STR_gene_cor(pairs, pt_str, norm_tumor, gene_cov, num_pc = 0):
    
    df_res = {
        "gene" : [],
        "str_id" : [],
        "coef" : [],
        "intercept" : [],
        "p_value" : []
    }

    for i in range(pairs.shape[0]):
        print(i)
        str_id, gene_id, gene_name = pairs.iloc[i,:]
        one_str = pt_str.loc[str_id,].to_frame()
        one_gene = norm_tumor.loc[gene_name,].to_frame()
        test_pair = pd.merge(one_str, one_gene, left_on = "patient", right_on = one_gene.index)
        test_pair.dropna(inplace=True)
        test_pair = test_pair.merge(gene_cov, left_on = "patient", right_on = gene_cov.index, how = "left")
        
        str_freq = test_pair[str_id].value_counts().to_frame()
        if str_freq.shape[0] < 3:
            continue
        else:
            numeric_vars =[str_id] + list(gene_cov.columns.values[:num_pc+1])
            formula = gene_name + ' ~ ' + ' + '.join(numeric_vars) + ' + C(gender) + C(ad_pop) + C(year)'
            model = smf.ols(formula, data = test_pair)
            result = model.fit()
            
            df_res["gene"].append(gene_id)
            df_res["str_id"].append(str_id)
            df_res["coef"].append(result.params[str_id])
            df_res["intercept"].append(result.params["Intercept"])
            df_res["p_value"].append(result.pvalues[str_id])
            
    df_res = pd.DataFrame(df_res)
    df_res["adj_p"] = multitest.multipletests(df_res["p_value"], method = "fdr_bh")[1]

    return df_res

def reg_out_cov(estr, pt_str, norm_tumor, gene_cov, num_pc = 0):
    
    resid_coef = []
    resid_inter = []
    for _, gene_str in estr.iterrows():
        gene_id = gene_str["gene"]
        str_id = gene_str["str_id"]
        gene_name = gene_id.replace("-", "_").replace(".", "_")
        one_str = pt_str.loc[str_id,].to_frame()
        one_gene = norm_tumor.loc[gene_name,].to_frame()    
        test_pair = pd.merge(one_str, one_gene, left_on = "patient", right_on = one_gene.index)
        test_pair.dropna(inplace=True)
        test_pair = test_pair.merge(gene_cov, left_on = "patient", right_on = gene_cov.index, how = "left")

        # regressed out covariates
        numeric_vars = list(gene_cov.columns.values[:num_pc+1]) #' + '.join(numeric_vars) +
        formula = gene_name + ' ~ ' +  '+ C(gender) + C(ad_pop) + C(year)'
        model = smf.ols(formula, data = test_pair).fit()

        X = sm.add_constant(test_pair[str_id])
        model_resid = sm.OLS(model.resid, X).fit()
        resid_coef.append(model_resid.params[str_id])
        resid_inter.append(model_resid.params["const"])
    estr["resid_coef"] = resid_coef
    estr["resid_inter"] = resid_inter
    
    return estr


if __name__ == "__main__":
    main()