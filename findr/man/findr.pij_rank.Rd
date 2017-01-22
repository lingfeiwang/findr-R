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
\name{findr.pij_rank}
\alias{findr.pij_rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation probability with expression data
}
\description{
Inference of pairwise correlation A- - -B probability between gene sets {A} and {B}, with expression data only, by converting log likelihood ratios to probabilities per A.
}
\usage{
findr.pij_rank(dt, dt2, nodiag = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dt}{
Input matrix of expression levels of A in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt,ns).
}
  \item{dt2}{
Input matrix of expression levels of B in matrix format. Element [i,j] is the expression level of gene i of sample j. The matrix has dimension (nt2,ns).
}
  \item{nodiag}{
When A and B are the same, log likelihood ratio between alternative and null hypotheses gives infinity. To avoid its contamination in the conversion from log likelihood ratios into probabilities, users need to arrange data accordingly, when {A} and {B} are the same or when {A} is a subset of {B}. The top submatrix of B's expression data must be identical with A, and nodiag must be set to TRUE. Otherwise, in the default configuration, {A} and {B} should not have any intersection and nodiag = FALSE.
}
}
\value{
\item{p}{Output matrix of inferred probability of correlation A- - -B versus no correlation A    B in matrix format. Element [i,j] is the probability of gene i being correlated with gene j. The matrix has dimension (nt,nt2).}
}
\note{
The methodology is to first calculate the log likelihood ratio (LLR) between each pair of A and B, and then convert the LLRs into probabilities for each A against all Bs.

nt: Number of genes for A.

nt2: Number of genes for B.

ns: Number of samples.
}

\seealso{
findr.pij_gassist

findr.pij_gassist_trad

findr.pijs_gassist
}
\examples{
library(findr)
data(geuvadis)
# Library initialization
a=findr.lib(loglv=12)
ans=findr.pij_rank(geuvadis$dt,geuvadis$dt2)
}
