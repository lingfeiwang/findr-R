\name{findr.pij_rank_a}
\alias{findr.pij_rank_a}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict pairwise regulation probability with expression data
}
\description{
Inference of pairwise correlation A- - -B probability between gene sets {A} and {B}, with expression data only, by converting log likelihood ratios to probabilities per A.
}
\usage{
findr.pij_rank_a(dt, dt2, nodiag = FALSE)
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
findr.pijs_gassist_a

findr.pijs_gassist_tot
}
\examples{
#TBA
}
