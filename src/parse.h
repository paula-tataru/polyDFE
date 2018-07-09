/*
 * polyDFE v2.0: predicting DFE and alpha from polymorphism data
 * Copyright (c) 2018  Paula Tataru and Marco A.P. Franco
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: paula@birc.au.dk
 */

#ifndef PARSE_H_
#define PARSE_H_

#define FILE_NOT_FOUND 1
#define PARAM_NOT_FOUND -1
#define LINE_SHORT -3
#define LINE_LONG 3

#include <gsl/gsl_vector_double.h>

#include "likelihood.h"

int parse_data(char *filename, Params *p, unsigned *n);
int parse(int what, char *filename, int id, void *pv, int initial_values);

#endif
