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
\docType{data}
\name{geuvadis}
\alias{geuvadis}
\title{Part of expression and best eQTL data from GEUVADIS}
\format{A list with the following variables:
\itemize{
  \item dgmi: haplotypes of the best cis-eQTLs of 10 miRNAs
  \item dmi: expression levels of 10 miRNAs
  \item dgt: haplotypes of the best cis-eQTLs of 1000 genes
  \item dt: expression levels of 1000 genes that have cis-eQTLs
  \item dt2: expression levels of all 3000 genes
  \item dc: continuous causal anchors for demonstration purposes, simulated from adding continuous noise to dgmi
}}
\usage{
data(geuvadis)
}
\source{\href{http://www.geuvadis.org/}{GEUVADIS Consortium}}
\description{The dataset is downsampled from the Geuvadis study, including expression levels of 10 miRNAs and 3000 genes for 360 European individuals. Among them, all miRNAs and 1000 genes have significant cis-eQTLs, whose haplotypes are also included.
}
\references{
Lappalainen, T. et al. Transcriptome and genome sequencing uncovers functional variation in humans. Nature 501, 506–511 (2013).
}
\keyword{datasets}
