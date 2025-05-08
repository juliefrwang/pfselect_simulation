
# copy from https://github.com/msesia/knockoff-filter/blob/master/R/knockoff/R/knockoff_filter.R
knockoff.threshold <- function(W, fdr=0.10, offset=0) {
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}
W <- results$W_statistic_matrix[1,]
knockoff.threshold(W)
