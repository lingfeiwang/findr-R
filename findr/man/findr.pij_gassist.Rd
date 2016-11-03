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
\name{findr.pij_gassist*}
\alias{findr.pij_gassist}
\alias{findr.pij_gassist_trad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation probability with expression and eQTL data
}
\description{
Inference of regulation E(A)->A->B probability with expression data of A,B and best eQTL data for A as E(A). For highest precision, the user should not pre-select secondary target genes (B), but instead include as many B's of the same type as possible, and pick the gene pairs of interest only after Findr's calculation. This function does not include p1.

findr.pij\_gassist performs our recommeneded combination of tests.

findr.pij\_gassist\_trad performs traditional causal inference test, but is not and is not intended as a reimplementation of Trigger package.
}
\usage{
findr.pij_gassist(dg, dt, dt2, na = NULL, nodiag = FALSE)
findr.pij_gassist_trad(dg, dt, dt2, na = NULL, nodiag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dg}{
Input matrix of best eQTL genotype data E(A), each row of which is the best eQTL of the corresponding row of t. Data is in matrix format. Element [i,j] is the genotype value of the best eQTL of gene i of sample j, and should be among values 0,1,...,na+1. The matrix has dimension (nt,ns).
}
  \item{dt}{
Input matrix of expression levels of A in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt,ns).
}
  \item{dt2}{
Input matrix of expression levels of B in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt2,ns)
}
  \item{na}{
Number of alleles for the species considered. This constrains every genotype data to be among 0,1,...,na+1. If unspecified (NULL), na is automatically determined as the maximum value of dg+1.
}
  \item{nodiag}{
When A and B are the same, log likelihood ratio between alternative and null hypotheses gives infinity. To avoid its contamination in the conversion from log likelihood ratios into probabilities, users need to arrange data accordingly, when {A} and {B} are the same or when {A} is a subset of {B}. The top submatrix of B's expression data must be identical with A, and nodiag must be set to TRUE. Otherwise, in the default configuration, {A} and {B} should not have any intersection and nodiag = FALSE.
}
}
\value{
Output matrix of recommended combination of multiple test results. Element [i,j] is the test probability of A->B for A = gene i and B = gene j. The matrix has dimension (nt,nt2).
}
\note{
The methodology is to first calculate the log likelihood ratio (LLR) for each test, and then convert the LLRs into probabilities, for each A against all Bs.

nt: Number of genes for A.

nt2: Number of genes for B.

ns: Number of samples.
}
\seealso{
findr.pij_rank

findr.pijs_gassist
}
\examples{
library(findr)
data(geuvadis)
# Library initialization
a=findr.lib(loglv=12)
ans=findr.pij_gassist(geuvadis$dg,geuvadis$dt,geuvadis$dt2)
}
