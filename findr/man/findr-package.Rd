% Copyright 2016, 2017 Lingfei Wang
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
Fast Inference of Networks from Directed Regulations
}
\description{
This is an R package of findr library. Findr (Fast Inference of Networks from Directed Regulations) is a C library for the inference of gene regulatory network structure. Findr performs causal inferences of pairwise gene regulations based on genotype and gene expression datasets or from gene expression dataset alone. Findr then bases on the causal inference results as prior edge information to construct a directed acyclic graph. A sparse Bayesian network can then be inferred with the R package lassopv to produce the p-value of every possible regulation. Methodologies can be found in publication (Wang and Michoel, 2016, Wang and Michoel, 2017).
}
\details{

%The DESCRIPTION file:
%\packageDESCRIPTION{findr}
%\packageIndices{findr}

This is the R package of findr. It contains findr C library, part of GNU Scientific Library (GSL), and the R interface. Findr is a C library intended for the inference of gene regulatory network structures. Findr performs causal inferences of pairwise gene regulations based on genotype and gene expression datasets or from gene expression dataset alone. Findr then bases on the causal inference results as prior edge information to construct a directed acyclic graph. A sparse Bayesian network can then be inferred with the R package lassopv to produce the p-value of every possible regulation. Methodologies can be found in publication (Wang and Michoel, 2016, Wang and Michoel, 2017).

For a complete manual, see: https://github.com/lingfeiwang/findr/blob/master/doc.pdf.

The R package of findr should first be initialized before using any functionality with the function:

findr.lib



The p-value of pairwise gene regulation can be inferred with any of the following functions. This can be more accurate than probabilities when the number of secondary targets is small ( several hundred or less ):

findr.pij_rank_pv

findr.pijs_gassist_pv

findr.pijs_cassist_pv



The probability of pairwise gene regulation can be inferred with any of the following functions:

findr.pij_rank

findr.pij_gassist

findr.pij_gassist_trad

findr.pijs_gassist

findr.pij_cassist

findr.pij_cassist_trad

findr.pijs_cassist



The reconstruction of a directed acyclic graph, based on the output of above pairwise function, can obtained from the following function:

findr.netr_one_greedy
}
\author{
Lingfei Wang <Lingfei.Wang.github@outlook.com>

Maintainer:  Lingfei Wang <Lingfei.Wang.github@outlook.com>
}
\references{
Lingfei Wang and Tom Michoel (2017) Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data. PLOS Computational Biology 13(8): e1005703. https://doi.org/10.1371/journal.pcbi.1005703

Lingfei Wang and Tom Michoel, Scalable causal gene network inference with prior node ordering and posterior FDR control (submitted). 2017.
}
\keyword{package}
\examples{
library(findr)
data(geuvadis)
# Library initialization
a=findr.lib(loglv=12)

# Predict pairwise correlation without genotype data
ans1=findr.pij_rank(geuvadis$dt,geuvadis$dt2)

# Predict pairwise regulation with genotype data with recommended test combination
ans2=findr.pij_gassist(geuvadis$dg,geuvadis$dt,geuvadis$dt2)

# Predict pairwise regulation with genotype data with all tests,
# only if User needs to combine test results on their own from findr output.
ans3=findr.pijs_gassist(geuvadis$dg,geuvadis$dt,geuvadis$dt2)

# Or similarly with continuous anchor instead of eQTLs
ans4=findr.pij_cassist(geuvadis$dc,geuvadis$dt,geuvadis$dt2)

# Construct maximal directed acyclic graph based on pairwise edge significance
ans5=findr.netr_one_greedy(ans2[1:dim(ans2)[1],1:dim(ans2)[1]])

# P-values can also be obtained, allowing more accurate results with limited secondary 
# targets, and allowing downstream FDR control with other packages e.g. qvalue or fdrtool
ans6=findr.pijs_gassist_pv(geuvadis$dg,geuvadis$dt,geuvadis$dt2)
}
