/*
 * polyDFE v2.0: predicting DFE and alpha from polymorphism data
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

#ifndef SIM_H_
#define SIM_H_

#include <gsl/gsl_rng.h>

#include "likelihood.h"

int sim_data(ParamsModel *pm, double *len, char *filename);

#endif
