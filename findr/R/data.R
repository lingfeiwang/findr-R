#' Part of expression and best eQTL data from GEUVADIS
#'
#' @docType data
#'
#' @usage data(geuvadis)
#'
#' @source \href{http://www.geuvadis.org/}{GEUVADIS Consortium}
#'
#' @format A list with 3 variables:
#' \itemize{
#'   \item dt: miRNA expression data for targets A as a matrix (10,360)
#'   \item dt2: gene expression data for targets B as a matrix (3000,360)
#'   \item dg: genotype data for best eQTLs of targets A as a matrix (10,360)
#' }
"geuvadis"
