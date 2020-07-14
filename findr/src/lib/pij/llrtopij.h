/* Copyright 2016-2018, 2020 Lingfei Wang
 * 
 * This file is part of Findr.
 * 
 * Findr is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Findr is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with Findr.  If not, see <http://www.gnu.org/licenses/>.
 */
/* This file contains the conversion from log likelihood ratio to probabilities
 *
 */

#ifndef _HEADER_LIB_PIJ_LLRTOPIJ_H_
#define _HEADER_LIB_PIJ_LLRTOPIJ_H_
#include "../base/config.h"
#include "../base/gsl/histogram.h"
#include "../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif


/* Use central histogram to estimate distribution probabilities of any point
 * within the histogram range. Linear intepolation is used.
 * Points outside histogram range gives boundary output
 * hc:	central histogram for estimation
 * d:	data (x coordinates of histogram) to be estimated their probabilities
 * ans:	output of estimated probabilities
 */
void pij_llrtopij_histogram_interpolate_linear(const gsl_histogram *hc,const VECTORF* d,VECTORF* ans);
 
/* Calculate buffer sizes for histogram conversion in pij_llrtopij_convert_histograms_buffed.
 * n:	Number of histogram bins. This must match pij_llrtopij_convert_histograms_buffed.
 * n1,
 * n2:	Sizes of two buffers for VECTORD.
 */
void pij_llrtopij_convert_histograms_get_buff_sizes(size_t n,size_t *n1,size_t *n2);

/* Allocate buffer for histogram conversion in pij_llrtopij_convert_histograms_buffed.
 * n:	Number of histogram bins. This must match pij_llrtopij_convert_histograms_buffed.
 * vb1,
 * vb2:	Output locations of allocated buffers.
 * Return:	0 on success.
 */
int pij_llrtopij_convert_histograms_make_buffs(size_t n,VECTORD** vb1,VECTORD** vb2);

/* Convert density histograms of null and real distribution into probability central
 * histogram with buffer provided. Both histograms must be distributions
 * (sum to unity and nonnegative).
 * hreal:	(n) Real density histogram to convert from. Also changed in calculation.
 * vnull:	(n) Null density histogram in vector format. Also changed in calculation.
 * hc:		(n+2) Central probability histogram as output.
 * vb1,
 * vb2:		Buffers needed for conversion. To allocate buffers, use
 			pij_llrtopij_convert_histograms_make_buffs.
 */
void pij_llrtopij_convert_histograms_buffed(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc,VECTORD* vb1,VECTORD* vb2);

/* Convert density histograms of null and real distribution into probability central histogram. Both histograms must be distributions (sum to unity and nonnegative).
 * hreal:	(n) Real density histogram to convert from. Also changed in calculation.
 * vnull:	(n) Null density histogram in vector format. Also changed in calculation.
 * hc:		(n+2) Central probability histogram as output.
 * Return:	0 if success.
 */
int pij_llrtopij_convert_histograms(gsl_histogram* hreal,VECTORD* vnull,gsl_histogram* hc);


/* Obtains the maximum of matrix, possibly ignoring diagonal elements.
 * Fails in the presence of NAN, and warns and updates at INFs.
 * d:		Matrix/Vector to get maximum, and update any INFs
 * nodiag:	Whether to ignore diagonal values when searching for maximum.
 * Return:	0 if NAN is found, or the non-INF maximum otherwise.
 */
FTYPE pij_llrtopij_llrmatmax(MATRIXF* d,char nodiag);
FTYPE pij_llrtopij_llrvecmax(VECTORF* d);


/* Convert LLR of real data to probabilities, when the distribution
 * of LLR of null distribution can be calculated analytically to follow
 * x=-0.5*log(1-z1/(z1+z2)), where z1~chi2(n1),z2~chi2(n2).
 * The conversion is performed for each gene A, i.e. per row of d and dconv.
 * This function is older than pij_llrtopij_convert_single so it is not parallel.
 * Make it parallel before using.
 * d:		[nrow,nx] The data to use for calculation of conversion rule from LLR to pij.
 * dconv:	[nrow,nd] The data of LLR to actually convert to pij. Can be same with d.
 * ans:		[nrow,nd] The output location of converted pij from dconv.
 * n1,
 * n2:		Parameters of null distribution.
 * nodiag:	Whether diagonal elements of d should be ignored when converting
 * 			to probabilities.
 * nodiagshift:	Diangonal column shift for nodiag==1.
 * 				For nodiagshift>0/<0, use upper/lower diagonal. 
 * Return:	0 on success.
 */
int pij_llrtopij_convert_single(const MATRIXF* d,const MATRIXF* dconv,MATRIXF* ans,size_t n1,size_t n2,char nodiag,long nodiagshift);

// Same with pij_llrtopij_convert_single, for d=dconv=ans. Saves memory.
int pij_llrtopij_convert_single_self(MATRIXF* d,size_t n1,size_t n2,char nodiag,long nodiagshift);




#ifdef __cplusplus
}
#endif
#endif
