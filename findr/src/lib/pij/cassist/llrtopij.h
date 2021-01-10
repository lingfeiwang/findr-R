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
/* This file contains the conversion from log likelihood ratio to probabilities
 *
 */

#ifndef _HEADER_LIB_PIJ_CASSIST_LLRTOPIJ_H_
#define _HEADER_LIB_PIJ_CASSIST_LLRTOPIJ_H_
#include "../../base/config.h"
#include "../../base/types.h"
#ifdef __cplusplus
extern "C"
{
#endif

/* Converts four LLRs into probabilities together.
 * Uses pij_cassit_llrtopij1_a to pij_cassit_llrtopij5_a.
 * See above functions for parameter definitions.
 * Return: 0 if all functions are successful.
 */
int pij_cassist_llrtopijs(VECTORF* p1,MATRIXF* p2,MATRIXF* p3,MATRIXF* p4,MATRIXF* p5,size_t ns,char nodiag);














#ifdef __cplusplus
}
#endif
#endif
