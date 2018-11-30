
#' fitness
#'
#' Determine the fitness of a model run
#' @param results data.frame containing the results of a model run
#' @param penalties List containing penalty coefficients
#' @export
#' v

fitness <- function (results = data.frame(),
                     penalties=list(
                     THETA=5,
                     ETA=5,
                     EPS=5,
                     cov.success=-10,
                     min.success=-20,
                     cov.warnings=0,
                     boundary=0,
                     rounding=0,
                     zero.grad=0,
                     final.zero.grad=0,
                     hessian.reset=0,
                     s.singular=0,
                     sig.digs=0,
                     condition.number=0)
){
  results$ofv+
    penalties$THETA           * results$NTHETA+
    penalties$ETA             * results$NETA+
    penalties$EPS             * results$NEPS+
    penalties$cov.success     * results$covariance_step_successful+
    penalties$min.success     * results$minimization_successful+
    penalties$cov.warnings    * results$covariance_step_warnings+
    penalties$boundary        * results$estimate_near_boundary+
    penalties$rounding        * results$rounding_errors+
    penalties$zero.grad       * results$zero_gradients+
    penalties$final.zero.grad * results$final_zero_gradients+
    penalties$hessian.reset   * results$hessian_reset+
    penalties$s.singular      * results$s_matrix_singular

}