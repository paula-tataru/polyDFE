/*
 * polyDFE v1.0: predicting DFE and alpha from polymorphism data
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

#include "likelihood.h"

#include <float.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_pow_int.h>
#include <gsl/gsl_sys.h>
#include <math.h>
#include <stdlib.h>

#include "localoptim.h"
#include "transform.h"

/****************************************************************************
 * Functions on structures: printing
 ****************************************************************************/
void fprintf_groups(ParamsModel pm, FILE *f, char *s)
{
    unsigned i;
    unsigned curr_group = 1;
    fprintf(f, "%s Using %d groups for the r parameters: [", s, pm.no_groups - 1);
    // skip r1
    for (i = 1; i < pm.no_r - 1; i++)
    {
        if (pm.inv_groups[i] != curr_group)
        {
            curr_group = pm.inv_groups[i];
            fprintf(f, "], [");
        }
        fprintf(f, "%d", i+1);
        if (pm.inv_groups[i + 1] == curr_group)
        {
            fprintf(f, ", ");
        }
    }
    // print last r
    i = pm.no_r - 1;
    if (pm.inv_groups[i] != curr_group)
    {
        curr_group = pm.inv_groups[i];
        fprintf(f, "], [");
    }
    fprintf(f, "%d] \n", i + 1);
}

void fprintf_params_model(ParamsModel pm, FILE *f, char *s)
{
    char *aux = malloc(sizeof(char) * 15);
    switch (pm.model)
    {
        case 1:
        {
            fprintf(f, "%s Model: A\n", s);
            break;
        }
        case 2:
        {
            fprintf(f, "%s Model: B\n", s);
            break;
        }
        case 3:
        {
            fprintf(f, "%s Model: C\n", s);
            break;
        }
        case 4:
        {
            fprintf(f, "%s Model: D\n", s);
            break;
        }
    }

    // print neut params
    fprintf(f, "%s", s);
    print_hearder_solution_neut(pm, TRUE, f);
    fprintf(f, "\n%s", s);
    print_solution_neut(pm, TRUE, NULL, aux, f);
    fprintf(f, "\n%s", s);

    // print sel params
    print_hearder_solution_sel(pm, TRUE, f);
    fprintf(f, "\n%s", s);
    print_solution_sel(pm, TRUE, NULL, aux, f);
    fprintf(f, "\n%s", s);

    // print demo params
    print_hearder_solution_demo(pm, TRUE, f);
    fprintf(f, "\n%s", s);
    print_solution_demo(pm, TRUE, NULL, aux, f);
    fprintf(f, "\n");

    free(aux);
}

void print_expec(double *e, int n)
{
    unsigned i = 0;
    printf("Expec: \t");
    for (i = 0; i < n; i++)
    {
        printf("%.6f\t", e[i]);
    }
    printf("\n");
}

/****************************************************************************
 * Functions on structures: initialization
 * (set pointers to NULL)
 ****************************************************************************/
void initialize_params_model(ParamsModel *pm)
{
    pm->n = 0;
    pm->model = 0;
    pm->no_groups = 0;
    pm->no_r = 0;

    pm->r = NULL;
    pm->inv_groups = NULL;

    pm->eps_an = 0;
    pm->lambda = 0;

    // parameters for mutation distribution
    pm->theta_bar = 0;
    pm->a = 1;

    // parameters for the DFE
    pm->no_sel = 0;
    pm->sel_params = NULL;

    pm->div_flag = TRUE;

    // parameters ranges
    pm->r_min = 0;
    pm->r_max = 10;
    pm->eps_an_min = 0;
    pm->eps_an_max = 0.5;
    pm->lambda_min = 0;
    pm->lambda_max = 0.5;
    pm->theta_bar_min = 0;
    pm->theta_bar_max = 1;
    pm->a_min = 0;
    pm->a_max = 500;
    pm->sel_min = NULL;
    pm->sel_max = NULL;

    // parameters flags
    pm->r_flag = TRUE;
    pm->eps_an_flag = TRUE;
    pm->lambda_flag = TRUE;
    pm->theta_bar_flag = TRUE;
    pm->a_flag = TRUE;
    pm->sel_flag = NULL;
    pm->sel_fixed = 1;

    // number of parameters to estimate
    pm->neut = 0;
    pm->sel = 0;

    pm->k = 0.01;

    pm->w = NULL;
    pm->w_div = NULL;
}

void initialize_params(Params *p)
{
    p->counts_neut = NULL;
    p->counts_sel = NULL;
    p->no_neut = 0;
    p->no_sel = 0;
    p->expec_neut = NULL;
    p->expec_sel = NULL;

    p->verbose = 0;

    p->lnL = DBL_MAX;
    p->grad = DBL_MAX;

    p->counter = 0;
    p->max_counter = 5000;

    p->pm = malloc(sizeof(ParamsModel));
    initialize_params_model(p->pm);
}

void initialize_selection_params(ParamsModel *pm, char *range)
{
    // according to the model, set the right flags to TRUE
    unsigned i = 0;
    if (pm->model < 4)
    {
        for (i = 0; i < pm->no_sel; i++)
        {
            pm->sel_flag[i] = TRUE;
        }

        if (range == NULL)
        {
            // initialize ranges too
            // S_d
            pm->sel_min[0] = -100000;
            pm->sel_max[0] = 0;
            // b
            pm->sel_min[1] = 0;
            pm->sel_max[1] = 10;
            if (pm->model == 1)
            {
                // S_max
                pm->sel_min[2] = 0;
                pm->sel_max[2] = 100;
            }
            else
            {
                // p_b
                pm->sel_min[2] = 0;
                pm->sel_max[2] = 0.5;
                // S_b
                pm->sel_min[3] = 0;
                pm->sel_max[3] = 100;
            }
        }
    } else
    {
        // model D - every second flag should be TRUE
        for (i = 0; i < pm->no_sel / 2; i++)
        {
            // set the flag of the probabilities to TRUE
            pm->sel_flag[2 * i + 1] = TRUE;
            // set the flag of the coefficients to TRUE
            pm->sel_flag[2 * i] = FALSE;
        }

        // set the flag of the first parameter that is estimated
        // to be FALSE, as this parameter should be 1 - sum(rest)
        pm->sel_fixed = 1;
        pm->sel_flag[pm->sel_fixed] = 2;

        if (range == NULL)
        {
            // initialize ranges too
            for (i = 0; i < pm->no_sel; i++)
            {
                pm->sel_min[i] = 0;
                pm->sel_max[i] = 1;
            }
        }
    }
}

/****************************************************************************
 * Functions on structures: allocation
 * (and initialization of the allocated memory)
 ****************************************************************************/
void set_to_zero(double **arr, int len)
{
    unsigned i = 0;
    for (i = 0; i < len; i++)
    {
        (*arr)[i] = 0;
    }
}

void allocate_counts(Counts *c, int n)
{
    c->sfs = malloc(sizeof(double) * n);
    set_to_zero(&c->sfs, n);
}

int allocate_grouping(ParamsModel *pm, double *groups)
{
    unsigned i = 0;
    // with r for divergence; without r for divergence should be pm->n - 1
    pm->no_r = pm->n;

//    // set divergence without r
//    pm->no_r = pm->n - 1;

    // if I do not use divergence data, then I am working without r for divergence
    if (pm->div_flag == FALSE)
    {
        pm->no_r = pm->n - 1;
    }

    // allocate the inverse mapping
    pm->inv_groups = malloc(sizeof(int) * pm->no_r);

    // start with default grouping
    pm->no_groups = pm->no_r;
    for (i = 0; i < pm->no_r; i++)
    {
        pm->inv_groups[i] = i;
    }

    // read grouping from groups, if present
    if (groups)
    {
        // check that the groups are sorted
        for (i = 1; i < groups[1]; i++)
        {
            if (groups[i + 2] <= groups[i + 1])
            {
                return (WRONG_RANGE);
            }
        }

        // + 1 for r[1] = 1
        // + 1 because the number in file is 1 too small
        pm->no_groups = (int) groups[1] + 2;

        // first entry in r is always set to 1
        // keep it in its own group
        pm->inv_groups[0] = 0;
        // second entry in r is always in group 1
        pm->inv_groups[1] = 1;
        // update rest of the grouping info
        int curr_group = 1;
        for (i = 2; i < pm->no_r; i++)
        {
            if (groups[curr_group + 1] == i)
            {
                curr_group++;
            }
            pm->inv_groups[i] = curr_group;
        }

        // print info about the groups used to file
        fprintf_groups(*pm, stdout, "----");
    }

    // allocate r
    pm->r = malloc(sizeof(double) * pm->no_groups);
    set_to_zero(&pm->r, pm->no_groups);

    return (EXIT_SUCCESS);
}

void allocate_selection_params(ParamsModel *pm)
{
    // I also need to allocate the parameters for the DFE
    pm->sel_params = malloc(sizeof(double) * pm->no_sel);
    pm->sel_min = malloc(sizeof(double) * pm->no_sel);
    pm->sel_max = malloc(sizeof(double) * pm->no_sel);
    pm->sel_flag = malloc(sizeof(double) * pm->no_sel);

    // allocate integration workspace
    pm->w = gsl_integration_workspace_alloc(5000);
    pm->w_div = gsl_integration_workspace_alloc(5000);
}

void allocate_params_expec(Params *p)
{
    // requires that I know n
    p->expec_neut = malloc(sizeof(double) * p->pm->n);
    p->expec_sel = malloc(sizeof(double) * p->pm->n);

    set_to_zero(&p->expec_neut, p->pm->n);
    set_to_zero(&p->expec_sel, p->pm->n);
}

void allocate_params(Params *p)
{
    unsigned i = 0;

    p->counts_neut = malloc(sizeof(Counts) * p->no_neut);
    p->counts_sel = malloc(sizeof(Counts) * p->no_sel);

    for (i = 0; i < p->no_neut; i++)
    {
        allocate_counts(&p->counts_neut[i], p->pm->n);
    }
    for (i = 0; i < p->no_sel; i++)
    {
        allocate_counts(&p->counts_sel[i], p->pm->n);
    }

    allocate_params_expec(p);

    allocate_selection_params(p->pm);
}

/****************************************************************************
 * Functions on structures: freeing
 ****************************************************************************/
void free_counts(Counts *c)
{
    if (c->sfs)
    {
        free(c->sfs);
    }
}

void free_params_model(ParamsModel *pm)
{
    if (pm->r)
    {
        free(pm->r);
    }
    if (pm->inv_groups)
    {
        free(pm->inv_groups);
    }
    if (pm->sel_params)
    {
        free(pm->sel_params);
    }
    if (pm->sel_min)
    {
        free(pm->sel_min);
    }
    if (pm->sel_max)
    {
        free(pm->sel_max);
    }
    if (pm->sel_flag)
    {
        free(pm->sel_flag);
    }
    gsl_integration_workspace_free(pm->w);
    gsl_integration_workspace_free(pm->w_div);
}

void free_params(Params *p)
{
    unsigned i;
    for (i = 0; i < p->no_neut; i++)
    {
        free_counts(&p->counts_neut[i]);
    }
    for (i = 0; i < p->no_sel; i++)
    {
        free_counts(&p->counts_sel[i]);
    }
    if (p->counts_neut)
    {
        free(p->counts_neut);
    }
    if (p->counts_sel)
    {
        free(p->counts_sel);
    }
    if (p->expec_neut)
    {
        free(p->expec_neut);
    }
    if (p->expec_sel)
    {
        free(p->expec_sel);
    }
    free_params_model(p->pm);
    if (p->pm)
    {
        free(p->pm);
    }
}

/****************************************************************************
 * Functions on structures: other
 ****************************************************************************/
void copy_params_model(ParamsModel *pm, ParamsModel source)
{
    unsigned i = 0;

    pm->model = source.model;
    pm->n = source.n;
    pm->no_groups = source.no_groups;

    for (i = 0; i < pm->no_groups; i++)
    {
        pm->r[i] = source.r[i];
    }
    for (i = 0; i < pm->no_r; i++)
    {
        pm->inv_groups[i] = source.inv_groups[i];
    }

    pm->eps_an = source.eps_an;
    pm->lambda = source.lambda;
    pm->theta_bar = source.theta_bar;
    pm->a = source.a;

    pm->no_sel = source.no_sel;
    for (i = 0; i < pm->no_sel; i++)
    {
        pm->sel_params[i] = source.sel_params[i];
    }

    pm->r_min = source.r_min;
    pm->r_max = source.r_max;
    pm->eps_an_min = source.eps_an_min;
    pm->eps_an_max = source.eps_an_max;
    pm->lambda_min = source.lambda_min;
    pm->lambda_max = source.lambda_max;
    pm->theta_bar_min = source.theta_bar_min;
    pm->theta_bar_max = source.theta_bar_max;
    pm->a_min = source.a_min;
    pm->a_max = source.a_max;
    for (i = 0; i < pm->no_sel; i++)
    {
        pm->sel_min[i] = source.sel_min[i];
        pm->sel_max[i] = source.sel_max[i];
    }

    pm->k = source.k;

    pm->r_flag = source.r_flag;
    pm->eps_an_flag = source.eps_an_flag;
    pm->lambda_flag = source.lambda_flag;
    pm->theta_bar_flag = source.theta_bar_flag;
    pm->a_flag = source.a_flag;
    for (i = 0; i < pm->no_sel; i++)
    {
        pm->sel_flag[i] = source.sel_flag[i];
    }
    pm->sel_fixed = source.sel_fixed;

    pm->neut = source.neut;
    pm->sel = source.sel;
}

void set_params_sel(ParamsModel *pm, int no_steps, int *it)
{
    unsigned i = 0;
    double diff = 0;
	// move away from boundaries with per
	double per = 0.2;
	// sum of probabilities for model 4
	double sum = 1;

	// update the parameters
	for (i = 0; i < pm->no_sel; i++)
	{
	    if (pm->sel_flag[i] != FALSE)
	    {
	        diff = pm->sel_max[i] - pm->sel_min[i];
	        pm->sel_params[i] = pm->sel_min[i] + diff * per
	                        + diff * (1 - 2 * per) * it[i] / no_steps;
	    }
	    if (i != pm->sel_fixed)
	    {
	        sum -= pm->sel_params[i];
	    }
	}
	if (pm->model == 4)
	{
	    pm->sel_params[pm->sel_fixed] = sum;
	}

	// if model A, B or C and the range for b covers < 1 and > 1
	// I want half of the tries for b to be less than 1
    // and the other half to be larger than 1
	// b is found on position 1
	i = 1;
	if (pm->model < 4)
	{
	    if (pm->sel_flag[i] == TRUE && pm->sel_min[i] < 1 && pm->sel_max[i] > 1)
	    {
	        int this_step = no_steps / 2;
	        if (it[i] < this_step)
	        {
	            diff = 1 - pm->sel_min[i];
	            pm->sel_params[i] = pm->sel_min[i] + diff * per
	                            + diff * (1 - 2 * per) * it[i] / this_step;
	        }
	        else
	        {
	            int this_it = it[i] - this_step;
	            diff = pm->sel_max[i] - 1;
	            pm->sel_params[i] = 1 + diff * per
	                            + diff * (1 - 2 * per) * this_it / this_step;
	        }
	    }
	}
}

void set_params_eps(ParamsModel *pm, int no_steps, int it)
{
    double diff = 0;
    // move away from boundaries with per
    double per = 0.2;
    if (pm->eps_an_flag == TRUE)
    {
        diff = pm->eps_an_max - pm->eps_an_min;
        pm->eps_an = pm->eps_an_min + diff * per
                        + diff * (1 - 2 * per) * it / no_steps;
    }
}

void count_params(ParamsModel *pm)
{
    // set the number of parameters that need to be estimated
    pm->neut = 0;
    pm->sel = 0;
    int i;

    if (pm->lambda_flag == TRUE)
    {
        pm->neut++;
    }
    if (pm->theta_bar_flag == TRUE)
    {
        pm->neut++;
    }
    if (pm->a_flag == TRUE)
    {
        pm->neut++;
    }
    if (pm->eps_an_flag == TRUE)
    {
        pm->neut++;
    }
    if (pm->r_flag == TRUE)
    {
        pm->neut += pm->no_groups - 1;
    }

    // parameters for DFE
    for (i = 0; i < pm->no_sel; i++)
    {
        if (pm->sel_flag[i] == TRUE)
        {
            pm->sel++;
        }
    }
}

/****************************************************************************
 * Neutral expectations
 ****************************************************************************/
double get_neut_expec(ParamsModel pm)
{
    // expectation for SFS
    if (pm.i < pm.n - 1)
    {
        double corr = pm.theta_bar * pm.r[pm.inv_groups[pm.i]] / (pm.i + 1);
        return (corr);
    }
    // expectation for divergence
    if (pm.div_flag == FALSE)
    {
        // I do not estimate lambda and r_n, assume = 1
        return ((1 + pm.theta_bar / pm.n));
    }
    // with r for divergence
    if (pm.no_r == pm.n)
    {
        return ((pm.lambda + pm.r[pm.inv_groups[pm.n - 1]] * pm.theta_bar / pm.n));
    }
    // without r for divergence
    return (pm.lambda + pm.theta_bar / pm.n);
}

void set_neut_expec(ParamsModel pm, double **expec_neut)
{
    for (pm.i = 0; pm.i < pm.n; pm.i++)
    {
        (*expec_neut)[pm.i] = get_neut_expec(pm);
    }
}

/****************************************************************************
 * Selected expectations
 ****************************************************************************/
double div_expec(double x, void *pv)
{
    ParamsDiv *pd = (ParamsDiv *) pv;
    double result = gsl_sf_pow_int(x, pd->n - 1)
                              * (1 - exp(-pd->S * (1 - x))) / (1 - x);
    return (result);
}

double taylor_approx(int n, double S)
{
    double term = 1;
    int k = 0;
    double result = 0;
    for (k = 0; k < 11; k++)
    {
        term *= (-S) / (n + k);
        result += term / (k + 1);
    }
    return (-result);
}

double approx_div_int(ParamsModel *pm, double S)
{
    double result = 0;
    double err = 0;
    double status = 0;

    ParamsDiv pd;
    pd.n = pm->n;
    pd.S = S;

    gsl_function F;
    F.function = &div_expec;
    F.params = &pd;

    status = gsl_integration_qag(&F, 0, 1, 0, 1e-5, 1000, 3, pm->w_div, &result,
                                 &err);

    if (status)
    {
        // if integration fails, use Taylor approximation
        // this happens when S is very negative!
        result = taylor_approx(pm->n, S);
    }

    return (result);
}

double get_sel_expec_fix_sel(ParamsModel *pm, double S)
{
    // when S is close to 0, use the neutral expectations
    if (fabs(S) <= 1e-10)
    {
        return (get_neut_expec(*pm));
    }

    // need to treat S > 0 and S < 0 separately
    // for numerical stability
    double eS = 0;
    double meS = 0;

    if (S < 0)
    {
        eS = exp(S);
        meS = eS - 1;
    }
    else
    {
        eS = exp(-S);
        meS = 1 - eS;
    }

    // expectation for SFS
    if (pm->i < pm->n - 1)
    {
        double hyper_corr = 1;

        // compute the hypergeometric functions
        // TODO if get nan's and weird behavior again, it could be from hyperg
        if (S < 0)
        {
            // compute the hypergeometric functions
            hyper_corr = eS - gsl_sf_hyperg_1F1_int(pm->i + 1, pm->n, S);
        }
        // eS != 0 is too restrictive, try meS != 1
        else if (meS != 1)
        {
            // if eS is 0, the final value is 1, and it is already set to that
            hyper_corr = 1 - eS * gsl_sf_hyperg_1F1_int(pm->i + 1, pm->n, S);
        }

        double dnum = meS * (pm->i + 1) * (pm->n - (pm->i + 1));
        double corr = pm->theta_bar * pm->r[pm->inv_groups[pm->i]] * pm->n
                        * hyper_corr / dnum;
        return (corr);
    }

    double div = 0;

    // expectation for divergence
    if (pm->div_flag == FALSE)
    {
        // I do not account for miss-attributed polymorphism
        div = S;
    }
    else
    {
        div = pm->lambda * S;

        // add miss-attributed polymorphism
        // with r for divergence
        if (pm->no_r == pm->n)
        {
            div += pm->r[pm->inv_groups[pm->n - 1]] * pm->theta_bar
                            * approx_div_int(pm, S);
        }
        else
        {
            // without r for divergence - assume 1
            div += pm->theta_bar * approx_div_int(pm, S);
        }
    }

    if (S < 0)
    {
        return (div * eS / meS);
    }
    return (div / meS);
}

double get_disp_ref_gamma(double S, void *pv)
{
    ParamsModel *pm = (ParamsModel *) pv;
    double s_bar = pm->sel_params[0];
    double b = pm->sel_params[1];
    // s_max is 0 for models B and C
    double s_max = 0;
    if (pm->model == 1)
    {
        s_max = pm->sel_params[2];
    }
    // undisplace and unreflect both S and s_bar
    double phi = gsl_ran_gamma_pdf(s_max - S, b, (s_max - s_bar) / b);

    if (gsl_isinf(phi) == 1)
    {
        phi = DBL_MAX;
    }

    return (phi);
}

double get_integrand_disp_ref_gamma(double S, void *pv)
{
    return (get_disp_ref_gamma(S, pv)
            * get_sel_expec_fix_sel((ParamsModel *) pv, S));
}

double get_integrand_exp(double S, void *pv)
{
    ParamsModel *pm = (ParamsModel *) pv;

    double s_b = pm->sel_params[3];
    double ex = gsl_ran_exponential_pdf(S, s_b);
    if (gsl_isinf(ex) == 1)
    {
        ex = DBL_MAX;
    }

    return (ex * get_sel_expec_fix_sel(pm, S));
}

int integrate_sel_expec(ParamsModel *pm, int no_expec, double **expec, double p,
                        double xmin, double xmax, void *func)
{
    int status = EXIT_SUCCESS;
    double res, err;
    double rel_err = 1e-5;
    double high_rel_err = rel_err;
    int key = 3;
    int high_key = key;

    gsl_function F;
    F.function = func;
    F.params = pm;

    // update p
    if (gsl_isnan(p))
    {
        p = 1;
    }

    // integrate each expectation separately
    for (pm->i = 0; pm->i < no_expec; pm->i++)
    {

        status = gsl_integration_qag(&F, xmin, xmax, 0, rel_err, pm->w->limit,
                                     key, pm->w, &res, &err);
        if (status)
        {
            // try with higher error
            high_rel_err = rel_err;
            while (status == GSL_EROUND && high_rel_err < 1e-1)
            {
                high_rel_err += 1e-5;
                status = gsl_integration_qag(&F, xmin, xmax, 0, high_rel_err,
                                             pm->w->limit, key, pm->w,
                                             &res, &err);
            }
            // try with higher keys
            high_key = key;
            while (status == GSL_EROUND && high_key < 6)
            {
                high_key++;
                status = gsl_integration_qag(&F, xmin, xmax, 0, high_rel_err,
                                             pm->w->limit, high_key, pm->w,
                                             &res, &err);
            }
            // try to integrate for singularities
            if (status == GSL_ESING)
            {
                status = gsl_integration_qags(&F, xmin, xmax, 0, high_rel_err,
                                              pm->w->limit, pm->w, &res, &err);
            }
            // try non-adaptive
            if (status)
            {
                size_t eval = 0;
                high_rel_err = rel_err;
                status = gsl_integration_qng(&F, xmin, xmax, 0, high_rel_err,
                                             &res, &err, &eval);
                // try with higher error
                while (status == GSL_EROUND && high_rel_err < 1e-1)
                {
                    high_rel_err += 1e-5;
                    status = gsl_integration_qng(&F, xmin, xmax, 0,
                                                 high_rel_err, &res, &err,
                                                 &eval);
                }
            }
            if (status)
            {
                return (EXIT_FAILURE);
            }
        }

        (*expec)[pm->i] += p * res;
    }

    return (EXIT_SUCCESS);
}

int set_sel_expec(ParamsModel *pm, double **expec, unsigned negative_only)
{
    unsigned i = 0;
    unsigned no_expec = pm->n;
    if (pm->div_flag == FALSE)
    {
        no_expec = pm->n - 1;
    }

    set_to_zero(expec, pm->n);

    if (pm->model < 4)
    {
        // requires numerical integration over selection coefficients
        int status = EXIT_SUCCESS;
        double xmin, xmax;
        double almost_zero = 1e-10;

        // s_max is 0 for models B and C
        double s_max = 0;
        if (pm->model == 1)
        {
            s_max = pm->sel_params[2];
        }
        // p_b is not defined for model A
        double p_b = GSL_NAN;
        double s_b = 0;
        if (pm->model != 1)
        {
            p_b = pm->sel_params[2];
            s_b = pm->sel_params[3];
        }

        /*****************************************************************
         * Negative selection
         *****************************************************************/
        xmax = 0;
        // calculate the lower limit of integration
        // use the values in i = 0, those should be the largest for
        // negative selection
        pm->i = 0;
        xmin = -10;
        double current_value = get_integrand_disp_ref_gamma(xmin, pm);
        xmin = xmin * 2;
        double new_value = get_integrand_disp_ref_gamma(xmin, pm);
        while ((new_value >= current_value) || (new_value > almost_zero))
        {
            // if new_value and current_value are both equal to 0
            // and xmin is larger than the mean (in the negative direction)
            // then I can stop
            if (new_value == 0 && current_value == 0 && xmin < pm->sel_params[0])
            {
                break;
            }
            current_value = new_value;
            xmin *= 2;
            new_value = get_integrand_disp_ref_gamma(xmin, pm);
        }

        status = integrate_sel_expec(pm, no_expec, expec, 1 - p_b, xmin,
                                     xmax, &get_integrand_disp_ref_gamma);

        if (status != EXIT_SUCCESS)
        {
            return (EXIT_FAILURE);
        }

        if (negative_only == FALSE)
        {
            /*****************************************************************
             * Positive selection
             *****************************************************************/
            switch (pm->model)
            {
                case 1:
                {
                    xmax = s_max;
                    xmin = 0;
                    status = integrate_sel_expec(pm, no_expec, expec, p_b,
                                                 xmin, xmax,
                                                 &get_integrand_disp_ref_gamma);
                    if (status != EXIT_SUCCESS)
                    {
                        return (EXIT_FAILURE);
                    }
                    break;
                }
                case 2:
                {
                    if (p_b <= 0)
                    {
                        break;
                    }
                    for (pm->i = 0; pm->i < no_expec; pm->i++)
                    {
                        (*expec)[pm->i] += p_b
                                            * get_sel_expec_fix_sel(pm, s_b);
                    }
                    break;
                }
                case 3:
                {
                    if (p_b <= 0)
                    {
                        break;
                    }
                    xmin = 0;
                    // calculate the upper limit of integration
                    // use the values in i = n - 1, those should be the largest for
                    // positive selection
                    // I could use i = n (divergence), but lambda reduces the value
                    pm->i = pm->n - 1;
                    xmax = 10;
                    while (get_integrand_exp(xmax, pm) > almost_zero)
                    {
                        xmax *= 2;
                    }

                    status = integrate_sel_expec(pm, no_expec, expec, p_b, xmin,
                                                 xmax, &get_integrand_exp);

                    if (status != EXIT_SUCCESS)
                    {
                        return (EXIT_FAILURE);
                    }
                    break;
                }
            }
        }
    }
    else
    {
        // discretized selection
        // I can, in theory, optimize this part of the code
        // as the expensive part of the expectation will remain constant
        // but this will involve quite a bit of recoding...
        // and it is rather fast anyways
        for (i = 0; i < pm->no_sel / 2; i++)
        {
            for (pm->i = 0; pm->i < no_expec; pm->i++)
            {
                if (negative_only == FALSE ||
                    (negative_only == TRUE && pm->sel_params[2*i] <= 0))
                (*expec)[pm->i] += pm->sel_params[2*i+1]
                          * get_sel_expec_fix_sel(pm, pm->sel_params[2*i]);
            }
        }
    }

    return (EXIT_SUCCESS);
}

double get_anc_expec(double eps_an, int i, int n, double *e)
{
    // combine the original calculated expectations
    // which do not contain the ancestral error
    // do this for entry i
    if (i < n - 1)
    {
        return ((1 - eps_an) * e[i] + eps_an * e[n - (i + 2)]);
    }
    return (e[i]);
}

void set_anc_expec(double eps_an, int n, double **e)
{
    // combine the original calculated expectations
    // which do not contain the ancestral error
    // do this for the full array e
    unsigned i = 0;
    double aux1 = 0;
    double aux2 = 0;

    for (i = 0; i < (n - 1)/2; i++)
    {
        aux1 = get_anc_expec(eps_an, i, n, *e);
        aux2 = get_anc_expec(eps_an, n - (i + 2), n, *e);
        (*e)[i] = aux1;
        (*e)[n - (i + 2)] = aux2;
    }
}

/****************************************************************************
 * Likelihood calculation for one fragment
 ****************************************************************************/
double get_log_neg_bin(int x, double mu, double a)
{
    // mu should never be 0
    if (mu == 0)
    {
        return (GSL_NAN);
    }

    // if a = -1, this is a sign I should use the Poisson distribution
    if (a == -1)
    {
        if (x == 0)
        {
            return (- mu);
        }

        return (-mu + x * log(mu) - gsl_sf_lngamma(x + 1));
//        return (log(gsl_ran_poisson_pdf(x, mu)));
    }

    // transform mu to probability
    if (x == 0)
    {
        return (a * (log(a) - log(a + mu)));
    }

    return (-log(x) - gsl_sf_lnbeta(x, a) + a * log(a) + x * log(mu)
                    - (a + x) * log(a + mu));
}

double get_lnL_one_frag(int div_flag, int n, Counts c, double *e, double a)
{
    double lnL = 0.0;
    unsigned i = 0;

    for (i = 0; i < n - 1; i++)
    {
        lnL += get_log_neg_bin(c.sfs[i], c.len_poly * e[i], a);

    }
    if (div_flag == TRUE)
    {
        lnL += get_log_neg_bin(c.sfs[n - 1], c.len_div * e[n - 1], a);
    }

    return (lnL);
}

/****************************************************************************
 * Full likelihood calculation for neutral/selected counts
 ****************************************************************************/
double get_lnL_aux(Params *p, int no, Counts *c, double *e)
{
    double lnL = 0.0;
    double curr_ln = 0.0;
    unsigned i;

    for (i = 0; i < no; i++)
    {
        curr_ln = get_lnL_one_frag(p->pm->div_flag, p->pm->n, c[i], e, p->pm->a);
        if (gsl_isnan(curr_ln) == 1)
        {
            break;
        }
        lnL += curr_ln;
    }

    if (i < no)
    {
        // finished too early - means I found NAN
        return GSL_NAN;
    }

    return (lnL);
}

void set_sel_lnL(Params *p)
{
    int status = EXIT_SUCCESS;

    status = set_sel_expec(p->pm, &p->expec_sel, FALSE);
    if (status == EXIT_FAILURE)
    {
        p->lnL_sel = GSL_NAN;
    }
    else
    {
        // add ancestral and contamination error
        set_anc_expec(p->pm->eps_an, p->pm->n, &p->expec_sel);
        p->lnL_sel = get_lnL_aux(p, p->no_sel, p->counts_sel, p->expec_sel);
    }
}

void set_neut_lnL(Params *p)
{
    // calculate expectations
    // it relies on having already calculated the selected expectations
    set_neut_expec(*p->pm, &p->expec_neut);
    // add ancestral error
    set_anc_expec(p->pm->eps_an, p->pm->n, &p->expec_neut);
    p->lnL_neut = get_lnL_aux(p, p->no_neut, p->counts_neut, p->expec_neut);
}

double get_lnL(const gsl_vector *x, void *pv)
{
    Params *p = (Params *) pv;

    // check the counter
    if (p->counter >= p->max_counter)
    {
        // reached maximum number of allowed likelihood evaluations
        // per optimization iteration
        // so GSL might be stuck - force it to stop by returning Inf
        return (DBL_MAX);
    }

    // update counter
    p->counter++;

    // update the content of p from x
    undo_transform(p->pm, x);

    // check the validity of probabilities for discretized DFE
    if (p->pm->model == 4
                    && (p->pm->sel_params[p->pm->sel_fixed] < 0
                    || p->pm->sel_params[p->pm->sel_fixed] > 1))
    {
        // return (GSL_NAN);
        return (DBL_MAX);
    }

    set_sel_lnL(p);
    if (gsl_isnan(p->lnL_sel) == 1)
    {
        // return (GSL_NAN);
        return (DBL_MAX);
    }

    set_neut_lnL(p);

    return (-(p->lnL_neut + p->lnL_sel));
}
