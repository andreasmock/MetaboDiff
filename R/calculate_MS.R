#' Calculate module significance
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @return adds module significance measures to metadata slot
#' @param group_factor name of column in colData for grouping
#' @examples
#' calculare_MS(met_example, group_factors = c("tumor_normal","random_gender"))
#' @export
calculate_MS = function(met, group_factors){
   for (i in 1:length(group_factors)){
        df = metadata(met)[[paste0("ttest_",group_factors[i])]]
        df$modules = metadata(met)$modules
        res_df = plyr::ddply(df,"modules",plyr::summarise,av_adj_pval=mean(adj_pval))
        metadata(met)[[paste0("MS_",group_factors[i])]] = res_df
    }
    met
}
