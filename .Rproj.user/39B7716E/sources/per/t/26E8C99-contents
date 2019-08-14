# List of Auxiliary Functions ---------------------------------------------
# (1) check_input : check symmetric matrix / dist object

# (1) check_input ---------------------------------------------------------
#' @keywords internal
#' @noRd
check_input <- function(x){
  if (inherits(x, "dist")){
    x = as.matrix(x)
  } else {
    if ((!is.matrix(x))||(nrow(x)!=ncol(x))||(!isSymmetric(x))){
      stop("* CHECK_INPUT : input should be a symmetric matrix.")
    }
  }
  if (any(abs(diag(x)) > .Machine$double.eps)){
    stop("* CHECK_INPUT : diagonal elements should all be zeros.")
  }
  return(x)
}
