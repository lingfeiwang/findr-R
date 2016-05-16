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
\name{findr.lib}
\alias{findr.lib}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Findr library initialization
}
\description{
This function loads findr library and initializes it with provided log level, initial random seed, and maximum thread count.
}
\usage{
findr.lib(loglv = 6, rs = 0, nth = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{loglv}{
Log level of library findr. Only messages whose levels are not greater than the designated log level will be printed. Level 0: critical, 1-3: errors, 4-6: warnings, 7-9: information, 10-12: debug. Default value is level 6.
}
  \item{rs}{
Initial random seed. Zero means to use current time.
}
  \item{nth}{
Maximum number of parrallel threads in calculation. For best performance, this should not exceed the number of CPU cores. Default value zero indicates using the number of cores automatically detected, which is not always accurate.
}
}
\examples{
# For more detailed log
library(findr)
findr.lib(loglv = 12)
# Single threaded
findr.lib(nth = 1)
}
