#' Perform comparative analysis using T-Test
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @param group_factors character vector of group factors
#' @return adds results from comparative analysis to metadata slot
#' @examples
#' diff_test(met_example, group_factors = c("tumor_normal","random_gender"))
#' @export
diff_test <- function(met, group_factors) {
    metadata(met) = vector("list",0)
    for (i in 1:length(group_factors)){
        df = genefilter::rowttests(assays(met)[["norm_imputed"]],
                       fac = as.factor(colData(met)[[group_factors[i]]]))
        df_ihw = IHW::as.data.frame(IHW::ihw(df$p.value,
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

