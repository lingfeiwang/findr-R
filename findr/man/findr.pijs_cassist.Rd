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
\name{findr.pijs_cassist}
\alias{findr.pijs_cassist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation probability with expression data and continuous causal anchors
}
\description{
Inference of regulation E(A)->A->B probability with expression data of A,B and continuous anchor data for A as E(A). This is similar with using eQTL but treats E(A) as a continuous variable. For highest precision, the user should not pre-select secondary target genes (B), but instead include as many B's of the same type as possible, and pick the gene pairs of interest only after Findr's calculation.
}
\usage{
findr.pijs_cassist(dc, dt, dt2, nodiag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dc}{
Input matrix of continuous anchor data E(A), each row of which is the anchor of the corresponding row of dt. Data is in matrix format. Element [i,j] is the value of the continuous anchor of gene i of sample j. The matrix has dimension (nt,ns).
}
  \item{dt}{
Input matrix of expression levels of A in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt,ns).
}
  \item{dt2}{
Input matrix of expression levels of B in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt2,ns)
}
  \item{nodiag}{
When A and B are the same, log likelihood ratio between alternative and null hypotheses gives infinity. To avoid its contamination in the conversion from log likelihood ratios into probabilities, users need to arrange data accordingly, when {A} and {B} are the same or when {A} is a subset of {B}. The top submatrix of B's expression data must be identical with A, and nodiag must be set to TRUE. Otherwise, in the default configuration, {A} and {B} should not have any intersection and nodiag = FALSE.
}
}
\value{
\item{p1}{Output vector of inferred probability of test 1, E(A)->A (alternative) versus E(A) A (null), in array format. Element [i] is the probability of best eQTL of gene i regulates gene i. The vector has dimension (nt). For nodiag=F, because the function expects significant eQTLs, p1 always return 1. For nodiag=T, uses diagonal elements of p2. Consider replacing p1 with your own (1-FPR) from eQTL discovery.}
\item{p2}{Output matrix of inferred probability of test 2, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A<-B (null), in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p3}{Output matrix of inferred probability of test 3, E(A)->A->B (null) versus E(A)->A- - -B with E(A)->B (alternative), in matrix format. Element [i,j] is the probability of null hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p4}{Output matrix of inferred probability of test 4, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A B, in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p5}{Output matrix of inferred probability of test 5, E(A)->A- - -B with E(A)->B (alternative) versus B<-E(A)->A, in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
}
\note{
The methodology is to first calculate the log likelihood ratio (LLR) for each test, and then convert the LLRs into probabilities, for each A against all Bs.

nt: Number of genes for A.

nt2: Number of genes for B.

ns: Number of samples.
}
\seealso{
findr.pij_rank

findr.pij_cassist

findr.pij_cassist_trad

findr.pijs_cassist_pv

findr.pijs_gassist

}
\examples{
library(findr)
data(geuvadis)
ans=findr.pijs_cassist(geuvadis$dc,geuvadis$dmi,geuvadis$dt2)
}
