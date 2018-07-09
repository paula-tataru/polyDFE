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

#ifndef BH_H_
#define BH_H_

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector_double.h>

#include "likelihood.h"
#include "localoptim.h"

#define MAX_ITER 29

struct params_basin_hop_s
{
    double max_same;
    double max_iter;
    double temp;

    /* Parameters for adaptive step taking */
    double step;        ///< current step size
    double accept_rate; ///< target acceptance rate
    int interval;       ///< how often to update the step size
    double factor;      ///< step size is changed by this factor upon each update
    int no_step;        ///< total number of steps
    int no_accept;      ///< total number of accepted steps

    /* Store best optimum found*/
    int size;           ///< number of estimated parameters
    gsl_vector *x;      ///< optimum parameters found
    double lnL;         ///< optimum joint likelihood found
    double criteria;    ///< gradient or size of optimum joint likelihood found
    double best_lnL;    ///< best likelihood ever encountered
};
typedef struct params_basin_hop_s ParamsBasinHop;

int check_lim(double x, double m, double M);
int check_lim_update(double *x, double *m, double *M, int flag, char *s);
void initialize_params_basin_hop(ParamsBasinHop *pb);
void set_default_max_iter(ParamsBasinHop *pb);
void allocate_params_basin_hop(ParamsBasinHop *pb, ParamsShare ps);
void free_params_basin_hop(ParamsBasinHop *pb);

int run_basin_hopping(ParamsBasinHop *pb,
                      const gsl_multimin_fdfminimizer_type *type,
                      ParamsOptim po, ParamsShare *ps);

#endif
