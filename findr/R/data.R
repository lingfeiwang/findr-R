# Copyright 2016 Lingfei Wang
# 
# This file is part of Findr.
# 
# Findr is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Findr is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with Findr.  If not, see <http://www.gnu.org/licenses/>.
# 
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
