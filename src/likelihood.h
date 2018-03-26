/*
 * polyDFE v1.11: predicting DFE and alpha from polymorphism data
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

#ifndef LN_H_
#define LN_H_

#include <gsl/gsl_vector_double.h>
#include <gsl/gsl_integration.h>
#include <stdio.h>

//WRONG_RANGE relies on no more than 100 000 lines in the input files
#define WRONG_RANGE 100000
#define WRONG_INIT 2
#define TRUE 1
#define FALSE 0

struct counts_s
{
    int len_poly;      ///< length of fragments
    int len_div;       ///< length of fragment for divergence counts
    double *sfs;       ///< observed SFS; last entry is the divergence
};
typedef struct counts_s Counts;

struct params_model_s
{
    unsigned n;         ///< sample size (number of chromosomes)
    unsigned model;     ///< model for distribution of selection: A(0), B(1), C(2)

    int *inv_groups;    ///< a map from i to r indexes (for grouping of r)
    double *r;          ///< demography nuisance parameters
    int no_groups;      ///< number of groups
    int no_r;           ///< number of r parameter (+1) used
                        ///< tells if I am using divergence r or not

    double eps_an;      ///< error in determining the ancestral allele
    double lambda;      ///< divergence parameter
    double theta_bar;   ///< mean of distribution of mutation rates
    double a;           ///< shape of distribution of mutation rates

    // parameters for the DFE
    // I store them in an array
    // model A: mean s_bar, shape b, maximum coefficient s_max
    // model B: mean deleterious s_d, shape b,
    //          probability to be beneficial p_b and coefficient beneficial s_b
    // model C: mean deleterious s_d, shape b,
    //          probability to be beneficial p_b and mean beneficial s_b
    // model D: coefficients s_i, followed by their probabilities, p_i
    //          s_i are fixed, only p_i are estimated
    int no_sel;
    double *sel_params;

    int div_flag;             ///< flag on whether to use divergence data or not

    /* Parameters ranges */
    double r_min;
    double r_max;
    double eps_an_min;
    double eps_an_max;
    double lambda_min;
    double lambda_max;
    double theta_bar_min;
    double theta_bar_max;
    double a_min;
    double a_max;
    double *sel_min;
    double *sel_max;

    double k;   ///< parameter for logistic transformation

    /* Flags for which parameters to estimate */
    int r_flag;
    int eps_an_flag;
    int lambda_flag;
    int theta_bar_flag;
    int a_flag;
    int *sel_flag;
    int sel_fixed;

    /* Number of parameters that need to be estimated - as given by flags */
    int neut;
    int sel;
    int i; // this keeps track of which sel SFS I integrate over

    int neut_ln; // do I optimize neutral likelihood?
    int sel_ln;  // do I optimize selected likelihood?

    gsl_integration_workspace *w;
    gsl_integration_workspace *w_div;

    int inital_estimation;
    int initial_values;
};
typedef struct params_model_s ParamsModel;

struct params_div_s
{
    int n;
    double S;
};
typedef struct params_div_s ParamsDiv;

struct params_s
{
    Counts *counts_neut;
    Counts *counts_sel;
    unsigned no_neut;    ///< number of neutrally evolving fragments
    unsigned no_sel;    ///< number of selected fragments
    double *expec_neut;
    double *expec_sel;

    ParamsModel *pm;

    double lnL;          ///< likelihood for the current set of parameters
    double lnL_neut;
    double lnL_sel;

    double grad;    ///< gradient of the best solution

    // kind of likelihood to optimize
    // 0: neutral, selected
    // 1: neutral, selected, joint
    // 2: joint
    int kind;

    int verbose;

    unsigned int counter;
    unsigned int max_counter;
};
typedef struct params_s Params;

void fprintf_params_model(ParamsModel pm, FILE *f, char *s);
void print_expec(double *e, int n);

void initialize_params_model(ParamsModel *pm);
void initialize_params(Params *p);
void initialize_selection_params(ParamsModel *pm, char *range);

void set_to_zero(double **arr, int n);
void allocate_selection_params(ParamsModel *pm);
int allocate_grouping(ParamsModel *pm, double *groups);
void allocate_params(Params *p);

void free_params_model(ParamsModel *pm);
void free_params(Params *p);

void copy_params_model(ParamsModel *pm, ParamsModel source);

void set_params_sel(ParamsModel *pm, int no_steps, int *it);
void set_params_eps(ParamsModel *pm, int no_steps, int it);

void count_params(ParamsModel *pm);

int set_sel_expec(ParamsModel *pm, double **expec, unsigned negative_only);
void set_neut_expec(ParamsModel pm, double **expec_neut);
void set_anc_expec(double eps_an, int n, double **expec);

void set_sel_lnL(Params *p);
double get_sel_lnL(const gsl_vector *x, void *pv);
void set_neut_lnL(Params *p);
double get_neut_lnL(const gsl_vector *x, void *pv);
double get_lnL(const gsl_vector *x, void *pv);

#endif
