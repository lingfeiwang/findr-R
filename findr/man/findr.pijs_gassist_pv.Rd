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
\name{findr.pijs_gassist_pv}
\alias{findr.pijs_gassist_pv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation p-values with expression and eQTL data
}
\description{
Inference of regulation E(A)->A->B p-values with expression data of A,B and best eQTL data for A as E(A). For highest precision, the user should not pre-select secondary target genes (B), but instead include as many B's of the same type as possible, and pick the gene pairs of interest only after Findr's calculation.
}
\usage{
findr.pijs_gassist_pv(dg, dt, dt2, na = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dg}{
Input matrix of best eQTL genotype data E(A), each row of which is the best eQTL of the corresponding row of t. Data is in matrix format. Element [i,j] is the genotype value of the best eQTL of gene i of sample j, and should be among values 0,1,...,na. The matrix has dimension (nt,ns).
}
  \item{dt}{
Input matrix of expression levels of A in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt,ns).
}
  \item{dt2}{
Input matrix of expression levels of B in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt2,ns)
}
  \item{na}{
Number of alleles for the species considered. This constrains every genotype data to be among 0,1,...,na. If unspecified (NULL), na is automatically determined as the maximum value of dg.
}
}
\value{
\item{p1}{Output vector of p-values of test 1, E(A)->A (alternative) versus E(A) A (null), in array format. Element [i] is the probability of best eQTL of gene i regulates gene i. The vector has dimension (nt).}
\item{p2}{Output matrix of p-values of test 2, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A<-B (null), in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p3}{Output matrix of p-values of test 3, E(A)->A->B (null) versus E(A)->A- - -B with E(A)->B (alternative), in matrix format. Element [i,j] is the probability of null hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p4}{Output matrix of p-values of test 4, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A B, in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p5}{Output matrix of p-values of test 5, E(A)->A- - -B with E(A)->B (alternative) versus B<-E(A)->A, in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
}
\note{
The methodology is to first calculate the log likelihood ratio (LLR) for each test, and then convert the LLRs into p-values.

nt: Number of genes for A.

nt2: Number of genes for B.

ns: Number of samples.
}
\seealso{
findr.pijs_gassist

findr.pijs_cassist_pv

findr.pij_rank_pv
}
\examples{
library(findr)
data(geuvadis)
ans=findr.pijs_gassist_pv(geuvadis$dg,geuvadis$dmi,geuvadis$dt2)
}
