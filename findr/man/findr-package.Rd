% Copyright 2016-2018, 2020 Lingfei Wang
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
This is an R package of findr library. Findr (Fast Inference of Networks from Directed Regulations) is a C library for the inference of gene regulatory network structure. Findr performs causal inferences of pairwise gene regulations based on genotype and gene expression datasets or from gene expression dataset alone. Findr then bases on the causal inference results as prior edge information to construct a directed acyclic graph. A sparse Bayesian network can then be inferred with the R package lassopv to produce the p-value of every possible regulation. Method details can be found in publication (Wang and Michoel, 2017, Wang and Michoel, 2017).
}
\details{

%The DESCRIPTION file:
%\packageDESCRIPTION{findr}
%\packageIndices{findr}

This is the R package of findr. It contains findr C library, part of GNU Scientific Library (GSL), and the R interface. Findr is a C library intended for the inference of gene regulatory network structures. Findr performs causal inferences of pairwise gene regulations based on genotype and gene expression datasets or from gene expression dataset alone. Findr then bases on the causal inference results as prior edge information to construct a directed acyclic graph. A sparse Bayesian network can then be inferred with the R package lassopv to produce the p-value of every possible regulation. Methodologies can be found in publication (Wang and Michoel, 2017, Wang and Michoel, 2017).

For a complete manual, see: https://github.com/lingfeiwang/findr/blob/master/doc.pdf.

The major function of Findr is to reconstruct causal gene networks using genome-transcriptome variation data. The core function is findr.pij_gassist. Below is a complete list of functions.

The probability of pairwise gene regulation can be inferred with any of the following functions:

findr.pij_rank

findr.pij_gassist

findr.pij_gassist_trad

findr.pijs_gassist

findr.pij_cassist

findr.pij_cassist_trad

findr.pijs_cassist



The p-value of pairwise gene regulation can be inferred with any of the following functions. This can be more accurate than probabilities when the number of secondary targets is small ( several hundred or less ):

findr.pij_rank_pv

findr.pijs_gassist_pv

findr.pijs_cassist_pv



The reconstruction of a directed acyclic graph, based on the output of above pairwise function, can obtained from the following function:

findr.netr_one_greedy



The R package of findr can be re-initialized for custom verbosity, thread count, and initial random seed with the function:

findr.lib


}
\author{
Lingfei Wang <Lingfei.Wang.github@outlook.com>

Maintainer:  Lingfei Wang <Lingfei.Wang.github@outlook.com>
}
\references{
Lingfei Wang and Tom Michoel (2017) Efficient and accurate causal inference with hidden confounders from genome-transcriptome variation data. PLOS Computational Biology 13(8): e1005703. https://doi.org/10.1371/journal.pcbi.1005703

Lingfei Wang and Tom Michoel, Scalable causal gene network inference with prior node ordering and posterior FDR control (in preparation). 2017.
}
\keyword{package}
\examples{
library(findr)
data(geuvadis)

# Perform the novel causal inference test between all pairs of 1000 genes, using cis-eQTLs as causal anchors
ans1=findr.pij_gassist(geuvadis$dgt,geuvadis$dt,geuvadis$dt,nodiag=TRUE)

# Then constructs maximal directed acyclic graph from pairwise causal inference test probability
ans2=findr.netr_one_greedy(ans1)

}
