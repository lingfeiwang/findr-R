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
\name{findr.netr_one_greedy}
\alias{findr.netr_one_greedy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reconstruct directed acyclic graphs with greedy algorithm
}
\description{
Reconstruction of a single directed acyclic graph from the given edge significance prior by adding the most significant edge one at a time, whilst avoiding cycles. As an example, the prior edge significant can be obtained from pairiwise inference from findr.pij_gassist or findr.pij_rank.
}
\usage{
findr.netr_one_greedy(dp,namax=NULL,nimax=NULL,nomax=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dp}{
Input matrix of significance levels of all edges. Element [i,j] is the significance level of edge i to j. The matrix has prior dimension (nt,nt).
}
  \item{namax}{
Maximum total number of edges in the reconstructed network. Defaults to unlimited, i.e. nt*(nt-1)/2.
}
  \item{nimax}{
Maximum number of incoming edges per node in the reconstructed network. Defaults to unlimited.
}
  \item{nomax}{
Maximum number of outgoing edges per node in the reconstructed network. Defaults to unlimited.
}
}
\value{
Output logical matrix of reconstructed network in matrix format. Element [i,j]= TRUE if an edge exists from node i to j on the reconstructed network, and FALSE otherwise. The matrix has dimension (nt,nt).
}
\note{
nt: Number of nodes.
}
\seealso{
findr.pij_gassist

findr.pij_rank
}
\examples{
library(findr)
data(geuvadis)
prior=findr.pij_gassist(geuvadis$dgt,geuvadis$dt,geuvadis$dt,nodiag=TRUE)
net=findr.netr_one_greedy(prior)
}
