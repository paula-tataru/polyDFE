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

#ifndef LOD_H_
#define LOD_H_

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector_double.h>
#include <stdio.h>

#include "likelihood.h"

#define SAME_STATE -3
#define SAME_STATE_CIR -4
#define RESTART -5
#define OUT_OF_TIME -6
#define FOUND_NAN -7
#define MAX_LK -8

struct params_optim_s
{
    double eps_abs;     ///< threshold for gradient norm
    double step_size;   ///< initial step size in line search
    double tol;         ///< accuracy of line minimization
    double max_iter;    ///< maximum allowed iterations
    double minutes;     ///< maximum allowed minutes for one full run
    int grad_descent;   ///< TRUE if I am using a gradient descent
};
typedef struct params_optim_s ParamsOptim;

void initialize_params_optim(ParamsOptim *po);

void print_header_solution_neut(ParamsModel pm, int i, int with_fixed, FILE *f,
                                int flag);
void print_header_solution_demo(ParamsModel pm, int i, int with_fixed, FILE *f,
                                int flag);
void print_header_solution_sel(ParamsModel pm, int i, int with_fixed, FILE *f,
                               int flag);
void print_header_solution(ParamsModel pm, int i, int with_fixed, int flag,
                           int neut_ln, int sel_ln);
void print_solution_neut(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag);
void print_solution_demo(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag);
void print_solution_sel(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag);
void print_solution(ParamsModel pm, int with_fixed, char *s, int flag,
                    int neut_ln, int sel_ln);

int optimize(const void *type, ParamsOptim po, ParamsShare *ps, char *s);

#endif
