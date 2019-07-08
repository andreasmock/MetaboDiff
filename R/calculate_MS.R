#' Calculate module significance
#'
#' @param met fully preprocessed MultiAssayExperiment object with slots "raw", "imputed", "norm" and "norm_imputed"
#' @return adds module significance measures to metadata slot
#' @param group_factor name of column in colData for grouping
#' @examples
#' calculate_MS(met_example, group_factors = c("tumor_normal","random_gender"))
#' @export
calculate_MS = function(met, group_factors){
    for (i in 1:length(group_factors)) {
        id = grep(group_factors[i],names(metadata(met)))[1]
        df = metadata(met)[[id]]
        df$modules = metadata(met)$modules
        res_df = plyr::ddply(df, "modules", plyr::summarise, av_pval=mean(pval),
                             av_adj_pval = mean(adj_pval))
        res_df$av_fold_change = rep(0,nrow(res_df))
        for(x in 1:nrow(df)){
            if(sum(df$modules==(x-1)&df$pval<0.05)>0){
            res_df$av_fold_change[x] = median(df[df$modules==(x-1)&df$pval<0.05,"fold_change"],na.rm=TRUE)
            }
        }
        metadata(met)[[paste0("MS_",group_factors[i])]] = res_df
    }
    met
}
