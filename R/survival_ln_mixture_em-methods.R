#' Visualizes the path of the EM algorithm
#' @param x A fitted `survival_ln_mixture_em` object.
#' 
#' @param ... Not used.
#'
#' @export
plot.survival_ln_mixture_em <- function(x, ...) {
  iter <- value <- var <- NULL
  rlang::check_dots_empty()
  
  if(!requireNamespace("plotly")) {
    tidyr::pivot_longer(x$em_iterations, 1:(ncol(x$em_iterations) - 1),
                        names_to = 'var') |> 
      dplyr::mutate(cat = ifelse(startsWith(var, 'beta'),
                                 'beta',
                                 ifelse(startsWith(var, 'phi'),
                                        'phi', 'eta'))) |> 
      ggplot2::ggplot() +
      ggplot2::geom_path(ggplot2::aes(x = iter, y = value, color = var)) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~cat, scales = 'free_y') +
      ggplot2::labs(y = 'Value',
                    x = 'Iteration') +
      ggplot2::guides(color = 'none')
  } else {
    (tidyr::pivot_longer(x$em_iterations, 1:(ncol(x$em_iterations) - 1),
                         names_to = 'var') |> 
       dplyr::mutate(cat = ifelse(startsWith(var, 'beta'),
                                  'beta',
                                  ifelse(startsWith(var, 'phi'),
                                         'phi', 'eta'))) |> 
       ggplot2::ggplot() +
       ggplot2::geom_path(ggplot2::aes(x = iter, y = value, color = var)) +
       ggplot2::theme_bw() +
       ggplot2::facet_wrap(~cat, scales = 'free_y') +
       ggplot2::labs(y = 'Value',
                     x = 'Iteration') +
       ggplot2::guides(color = 'none')) |> 
      plotly::ggplotly()
  }
}
