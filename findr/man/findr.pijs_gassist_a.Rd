\name{findr.pijs_gassist_*}
\alias{findr.pijs_gassist_a}
\alias{findr.pijs_gassist_tot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation probability with expression and eQTL data
}
\description{
Inference of regulation E(A)->A->B probability with expression data of A,B and best eQTL data for A as E(A).
}
\usage{
findr.pijs_gassist_a(dg, dt, dt2, na = NULL, nodiag = FALSE)
findr.pijs_gassist_tot(dg, dt, dt2, na = NULL, nodiag = FALSE)
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
\item{p1}{Output vector of inferred probability of test 1, E(A)->A (alternative) versus E(A) A (null), in array format. Element [i] is the probability of best eQTL of gene i regulates gene i. The vector has dimension (nt). Since the library expects significant eQTLs only, all elements of this variable are set to one artificially.}
\item{p2b}{Output matrix of inferred probability of test 2 bold, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A B, in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p2c}{Output matrix of inferred probability of test 2 conservative, E(A)->A- - -B with E(A)->B (alternative) versus E(A)->A<-B (null), in matrix format. Element [i,j] is the probability of alternative hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
\item{p3}{Output matrix of inferred probability of test 3, E(A)->A->B (null) versus E(A)->A- - -B with E(A)->B (alternative), in matrix format. Element [i,j] is the probability of null hypothesis for A = gene i and B = gene j. The matrix has dimension (nt,nt2).}
}
\note{
The methodology is to first calculate the log likelihood ratio (LLR) for each test, and then convert the LLRs into probabilities, either for each A against all Bs in findr.pijs_gassist_a, or altogether in findr.pijs_gassist_tot.

nt: Number of genes for A.

nt2: Number of genes for B.

ns: Number of samples.
}
\seealso{
findr.pij_rank_a
}
\examples{
#TBA
}
