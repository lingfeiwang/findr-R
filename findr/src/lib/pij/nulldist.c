/* Copyright 2016-2018, 2020, 2021 Lingfei Wang
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
#include "../base/config.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "../base/gsl/blas.h"
#include "../base/gsl/math.h"
#include "../base/logger.h"
#include "../base/macros.h"
#include "../base/histogram.h"
#include "../base/data_process.h"
#include "nulldist.h"
#pragma GCC diagnostic ignored "-Wunused-parameter"



/*************************************************************
 * Generic functions for any step
 *************************************************************/

void pij_nulldist_pdfs(const VECTORD* loc,VECTORD* ans,const void* param)
{
	const struct pij_nulldist_pdfs_param *p=param;
	size_t	nd=loc->size;
	size_t	i;

	//Part 1: (1-exp(-2*x))^((n1-2)/2)
	for(i=0;i<nd;i++)
		VECTORDF(set)(ans,i,log(-math_sf_expminusone(-2*VECTORDF(get)(loc,i))));
	VECTORDF(scale)(ans,(double)p->n1/2-1);
	//Part 2: exp(-n2*x)
	gsl_blas_daxpy(-(double)p->n2,loc,ans);
	//Part 3: 2*Gamma((n1+n2)/2)/(Gamma(n1/2)*Gamma(n2/2)
	VECTORDF(add_constant)(ans,M_LN2+math_sf_lngammahalf(p->n1+p->n2)-math_sf_lngammahalf(p->n1)-math_sf_lngammahalf(p->n2));
	//Final: all log
	for(i=0;i<nd;i++)
		VECTORDF(set)(ans,i,exp(VECTORDF(get)(ans,i)));
}

void pij_nulldist_calcpdf_buffed(long n1c,size_t n1d,long n2c,size_t n2d,const VECTORD* loc,MATRIXD* ans,VECTORD* vb2)
{
	size_t	nd=loc->size;
	size_t	i,j;
	
	assert(n1d&&n2d);
	assert(nd&&(ans->size2==nd)&&(vb2->size==nd));
	//Calculate vb2=log(1-exp(-2x))
	VECTORDF(memcpy)(vb2,loc);
	VECTORDF(scale)(vb2,-2);
	for(i=0;i<nd;i++)
		VECTORDF(set)(vb2,i,log(-math_sf_expminusone(VECTORDF(get)(vb2,i))));
	
	//Calculate ans[0] without nv-dependent coefficients
	//ans[0]=-n2d*x+(n1d/2-1)*log(1-exp(-2x))
	{
		VECTORDF(view) vv=MATRIXDF(row)(ans,0);
		VECTORDF(memcpy)(&vv.vector,loc);
		VECTORDF(scale)(&vv.vector,-(double)n2d);
		gsl_blas_daxpy((double)n1d/2.-1,vb2,&vv.vector);
	}

	//Calculate vb2=0.5*n1c*log(1-exp(-2x))+x*n2c
	VECTORDF(scale)(vb2,0.5*(double)n1c);
	gsl_blas_daxpy((double)n2c,loc,vb2);

	//Calculate ans[i] without nv-dependent coefficients
	for(i=1;i<ans->size1;i++)
	{
		VECTORDF(const_view) vv1=MATRIXDF(const_row)(ans,i-1);
		VECTORDF(view) vv2=MATRIXDF(row)(ans,i);
		VECTORDF(memcpy)(&vv2.vector,&vv1.vector);
		VECTORDF(add)(&vv2.vector,vb2);
	}

	//Include nv-dependent coefficients
	for(i=0;i<ans->size1;i++)
	{
		VECTORDF(view) vv=MATRIXDF(row)(ans,i);
		VECTORDF(add_constant)(&vv.vector,(FTYPE)(M_LN2+math_sf_lngammahalf((size_t)((long)i*(n1c-n2c)+(long)(n1d+n2d)))-math_sf_lngammahalf((size_t)((long)i*n1c+(long)n1d))-math_sf_lngammahalf((size_t)(-(long)i*n2c+(long)n2d))));
	}
	//Convert log pdf to pdf
	for(i=0;i<ans->size1;i++)
		for(j=0;j<ans->size2;j++)
			MATRIXDF(set)(ans,i,j,exp(MATRIXDF(get)(ans,i,j)));
}

static int pij_nulldist_calcpdf(long n1c,size_t n1d,long n2c,size_t n2d,const VECTORD* loc,MATRIXD* ans)
{
#define	CLEANUP	AUTOFREEVEC(vb)
		AUTOALLOCVECD(vb,loc->size,30000)
		if(!vb)
			ERRRET("Not enough memory.")
		pij_nulldist_calcpdf_buffed(n1c,n1d,n2c,n2d,loc,ans,vb);
		CLEANUP
		return 0;
#undef	CLEANUP
}

int pij_nulldist_hist_pdf(const double* restrict range,size_t nbin,double* restrict hist,size_t n1,size_t n2,size_t n)
{
#define CLEANUP	CLEANVECD(loc)CLEANVECD(val)
	VECTORD	*loc,*val;
	VECTORDF(view)	vvh=VECTORDF(view_array)(hist,nbin);
	MATRIXDF(view)	mvv;
	size_t		i;
	size_t		nsp;
	
	assert(n&&(n<10));
	nsp=(size_t)1<<(n-1);
	loc=VECTORDF(alloc)(nbin*nsp);
	val=VECTORDF(alloc)(nbin*nsp);
	if(!(loc&&val))
		ERRRET("Not enough memory.")
	
	//Construct bin ranges
	{
		VECTORDF(const_view) vvc=VECTORDF(const_view_array)(range,nbin+1);
		histogram_finer_central(&vvc.vector,loc,nsp);
	}
	
	mvv=MATRIXDF(view_vector)(val,1,val->size);
	//Calculate bin values
	if(pij_nulldist_calcpdf(0,n1,0,n2,loc,&mvv.matrix))
	{
		CLEANUP
		return 1;
	}
	
	//Shrink to output
	VECTORDF(set_zero)(&vvh.vector);
	for(i=0;i<nsp;i++)
	{
		VECTORDF(const_view)	vvc=VECTORDF(const_subvector_with_stride)(val,i,nsp,nbin);
		VECTORDF(add)(&vvh.vector,&vvc.vector);
	}
	//Scale to 1
	VECTORDF(scale)(&vvh.vector,1./(double)nsp);
	
	CLEANUP
	return 0;
#undef	CLEANUP
}






































