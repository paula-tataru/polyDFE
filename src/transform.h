/*
 * polyDFE v1.11: predicting DFE and alpha from polymorphism data
 * Copyright (c) 2018  Paula Tataru
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

#ifndef TRANS_H_
#define TRANS_H_

#include <gsl/gsl_vector_double.h>

#include "likelihood.h"

void transform(gsl_vector **x, ParamsModel pm);
void undo_transform(ParamsModel *pm, const gsl_vector *x);

void random_step(ParamsModel *pm, double step);

#endif
