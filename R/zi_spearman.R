#######################################################################################
#' scHOT implemented Pimentel's zero-inflated (ZI) correlation
#' in their R package, providing implementations
#' of the ZI spearman and tau rank correlations.
#' scMACS implements a slight modification of
#' scHOT's zero-inflated correlation measure, by
#' returning NAs in the cases where the correlation
#' is undefined, and modifying it to use a C-backed correlation program.
#' Both references are provided below and are
#' referenced in all documentations to indicate
#' their work in implementing these methods.
#'
#'
#' Description:
#' the weightedZISpearman function calculates weighted rho\*,
#' where rho\* is described in Pimentel et al (2009).
#' This association measure is defined for zero-inflated,
#' non-negative random variables.
#'
#' @title weightedZISpearman
#' @param w weight vector, values should be between 0 and 1
#' @param x x and y are non-negative data vectors
#' @param y x and y are non-negative data vectors
#' @param ZI boolean flag. When set to false, this will skip zero-inflation and just calculate the normal spearman
#' @return \code{numeric} weighted rho* association value between x and y
#'
#'
#' @references scHOT
#' Ghazanfar, S., Lin, Y., Su, X. et al. Investigating higher-order interactions in
#'   single-cell data with scHOT. Nat Methods 17, 799–806 (2020).
#'   https://doi.org/10.1038/s41592-020-0885-x
#'
#' Zero-Inflated Correlation
#' Pimentel, Ronald Silva, "Kendall's Tau and Spearman's Rho for
#'   Zero-Inflated Data" (2009). Dissertations. 721.
#'   https://scholarworks.wmich.edu/dissertations/721
#'
#' @internal
#' @noRd

weightedZISpearman <- function(x, y, w = 1, ZI = TRUE) {

  # needs the original values, not the ranks

  if (any(x < 0 | y < 0)) {
    stop("x and/or y values have negative values")
  }
  if (length(x) != length(y)) {
    stop("x and y should have the same length")
  }
  if (length(w) == 1) {
    w <- rep(w, length(x))
  }

  if(!ZI){

    spearmanCorr <- wCorr::weightedCorr(x = x, y = y, weights = w, method = "Spearman")
    return(spearmanCorr)
    
  }

  posx <- x > 0
  posy <- y > 0
  pospos <- posx & posy

  p_11 <- sum(w * pospos) / sum(w)
  p_00 <- sum(w * (!posx & !posy)) / sum(w)
  p_01 <- sum(w * (!posx & posy)) / sum(w)
  p_10 <- sum(w * (posx & !posy)) / sum(w)


  if (any(pospos) & p_11 > 0) {
    rho_11 <- wCorr::weightedCorr(x = x[pospos], y = y[pospos], weights = w[pospos], method = "Spearman")
  } else {
    print("Zero inflated Spearman correlation is undefined,
          returning NA")
    rho <- NA
    return(rho)
  }

  rho_star <- p_11 * (p_01 + p_11) * (p_10 + p_11) * rho_11 +
    3 * (p_00 * p_11 - p_10 * p_01)

  if (is.na(rho_star)) {
    print("Zero inflated Spearman correlation is undefined,
          returning NA")
    rho <- NA
    return(rho)
  }



  return(rho_star)
}


# ##' scHOT functions, supersceded by used of wCorr's weighted correlation function, which backends to C for faster calculations.
# weightedSpearman = function(x, y, w = 1) {
#
#   if (length(x) != length(y)) {
#     stop("x and y should have the same length")
#   }
#   if (length(w) == 1) {
#     w <- rep(w, length(x))
#   }
#
#   keep = w > 0
#
#   xr = rank(x[keep])
#   yr = rank(y[keep])
#   return(weightedPearson(x = xr, y = yr, w = w[keep]))
# }
#
# weightedPearson = function(x, y, w = 1) {
#
#   if (length(x) != length(y)) stop("data must be the same length")
#
#   if (length(w) == 1) {
#     w <- rep(w, length(x))
#   }
#
#   nw = sum(w)
#   wssx = nw * sum(w * (x^2)) - sum(w * x)^2
#   wssy = nw * sum(w * (y^2)) - sum(w * y)^2
#   wssxy = nw * sum(w * x * y) - sum(w * x) * sum(w * y)
#   wcor = wssxy/sqrt(wssx * wssy)
#   return(wcor)
# }
