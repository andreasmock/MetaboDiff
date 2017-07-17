#' Perform comparative analysis using T-Test
#'
#' @param met MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factors vector of group factors
#' @return MultiAssayExperiment object with additional metadata slot containing the results of the comparative analysis
#' @examples
#' diff_test(met_example, group_factors = c("tumor_normal","random_gender"))
diff_test <- function(met, group_factors) {
    metadata(met) = vector("list",0)
    for (i in 1:length(group_factors)){
        df = rowttests(assays(met)[["norm_imputed"]],
                       fac = as.factor(colData(met)[[group_factors[i]]]))
        df_ihw = as.data.frame(ihw(df$p.value,
                                   as.numeric(apply(assays(met)[["norm_imputed"]],1,var)),
                                   alpha = 0.05,
                                   nbins = 20))
        res_df = data.frame(pval=df_ihw$pvalue,
                            adj_pval=df_ihw$adj_pvalue,
                            fold_change=df$dm,
                            var=df_ihw$covariate)
        metadata(met)[[paste0("ttest_",group_factors[i])]] = res_df
    }
    met
}

