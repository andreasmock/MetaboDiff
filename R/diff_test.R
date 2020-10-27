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
        if(length(levels(as.factor(colData(met)[[group_factors[i]]])))>2){
        coeff = sapply(1:nrow(assays(met)[["norm_imputed"]]),
                     function(x) aov(assays(met)[["norm_imputed"]][x,]~as.factor(colData(met)[[group_factors[i]]]))$coefficient)
        xlevels = unlist(aov(assays(met)[["norm_imputed"]][1,]~as.factor(colData(met)[[group_factors[i]]]))$xlevels)
        res = sapply(1:nrow(assays(met)[["norm_imputed"]]),
                     function(x) summary(aov(assays(met)[["norm_imputed"]][x,]~as.factor(colData(met)[[group_factors[i]]]))))
        res_df = data.frame(pval=as.vector(sapply(sapply(res,"[",i=5),"[",i=1)),
                            adj_pval=p.adjust(as.vector(sapply(sapply(res,"[",i=5),"[",i=1)),method = "fdr"),
                            dm=coeff[2,])
        metadata(met)[[paste0("anova_",group_factors[i],"_",paste(xlevels,collapse = "_vs_"))]] = res_df
        } else {
        xlev = levels(as.factor(colData(met)[[group_factors[i]]]))
        df = genefilter::rowttests(assays(met)[["norm_imputed"]],
                                   fac = as.factor(colData(met)[[group_factors[i]]]))
        res_df = data.frame(metabolite=rownames(df),
                            pval=df$p.value,
                            adj_pval=p.adjust(df$p.value,method="fdr"),
                            dm=df$dm)
        metadata(met)[[paste0("ttest_",group_factors[i],"_",xlev[2],"_vs_",xlev[1])]] = res_df
        }
    }

    met
}

