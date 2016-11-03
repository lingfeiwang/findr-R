% Copyright 2016 Lingfei Wang
% 
% This file is part of Findr.
% 
% Findr is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% Findr is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with Findr.  If not, see <http://www.gnu.org/licenses/>.
% 
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

This is the R package of findr. It contains findr C library, part of GNU Scientific Library (GSL), and the R interface. Findr is a C library intended for the inference of gene regulatory network structures. Currently findr only contains the inference of pairwise regulation, with or without the assistance of genotype data. Functions to reconstruct gene regulation networks will be added in future releases. Findr has been shown to be more accurate but hundreds to thousands of times faster than existing methods (Wang and Michoel, 2016), from its methodology and implementation advances.

For a complete manual, see: https://github.com/lingfeiwang/findr/blob/master/doc.pdf.

The R package of findr should first be initialized before using any functionality with the function:

findr.lib



The probability of pairwise gene regulation can be inferred with any of the following functions:

findr.pij_rank

findr.pij_gassist

findr.pij_gassist_trad

findr.pijs_gassist
}
\author{
\packageAuthor{findr}

Maintainer: \packageMaintainer{findr}
}
\references{
Lingfei Wang and Tom Michoel. (In preparation). 2016.
}
\keyword{ package }
\examples{
library(findr)
data(geuvadis)
# Library initialization
a=findr.lib(loglv=12)
# Predict pairwise correlation without genotype data
ans1=findr.pij_rank(geuvadis$dt,geuvadis$dt2)
# Predict pairwise regulation with genotype data with recommended test combination
ans2=findr.pij_gassist(geuvadis$dg,geuvadis$dt,geuvadis$dt2)
# Predict pairwise regulation with genotype data with all tests. User needs to combine test results on their own from findr output.
ans3=findr.pijs_gassist(geuvadis$dg,geuvadis$dt,geuvadis$dt2)
}
