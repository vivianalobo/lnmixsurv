#' @export
burnin_and_thin <- function(post, burnin, thin) {
  len <- dim(post[[1]])[1]
  ret <- lapply(
    post,
    function(x) posterior::subset_draws(x, iteration = seq(from = burnin + 1, len, by = thin))
  )
  return(ret)
}
