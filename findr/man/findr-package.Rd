\name{findr-package}
\alias{findr-package}
\alias{findr}
\docType{package}
\title{
\packageTitle{findr}
}
\description{
\packageDescription{findr}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{findr}
\packageIndices{findr}

This is the R package of findr. It contains findr C library, part of GNU Scientific Library (GSL), and the R interface. Due to the lack of functionality in Rtools, this package can only be built on non-Windows platforms. Findr is a C library intended for the inference of gene regulatory network structures. Currently findr only contains the inference of pairwise regulation, with or without the assistance of genotype data. Functions to reconstruct gene regulation networks will be added in future releases. Findr has been shown to be more accurate but hundreds to thousands of times faster than existing methods (Wang and Michoel, 2016), from its methodology and implementation advances.

For a complete manual, see: https://github.com/lingfeiwang/findr/blob/master/doc.pdf.

The R package of findr should first be initialized before using any functionality with the function:

findr.lib



The probability of pairwise gene regulation can be inferred with any of the following functions:

findr.pij_rank_a

findr.pijs_gassist_a

findr.pijs_gassist_tot
}
\author{
\packageAuthor{findr}

Maintainer: \packageMaintainer{findr}
}
\references{
Lingfei Wang and Tom Michoel. Improved microRNA target predictions from genetics of gene expression data with Findr (submitted). 2016.
}
\keyword{ package }
\examples{
#TBA
}
