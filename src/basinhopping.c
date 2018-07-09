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

/* The basin hopping algorithm implemented here is a reimplementation of
 * python's scipy basinhopping
 */

#include "basinhopping.h"

#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sys.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "transform.h"

/****************************************************************************
 * Functions for estimating automated initial point
 ****************************************************************************/
int check_lim(double x, double m, double M)
{
    // verify that m <= x <= M
    if (x < m || x > M)
    {
        return (WRONG_RANGE);
    }
    return (EXIT_SUCCESS);
}

int check_lim_update(double *x, double *m, double *M, int flag, char *s)
{
    double aux = 0;

	if (flag == FALSE)
    {
        return (EXIT_SUCCESS);
    }

    // verify that m <= x <= M
    if (*x < *m || *x > *M)
    {
        printf("----      The initial value %g of %s did not fit the given [%g, %g] range.\n",
               *x, s, *m, *M);

        // update limits so that x is within the limits
        if (*x < *m)
        {
            // make m smaller
            aux = (*M - *m) * 0.1;
        	*m = *x - aux;
        	// if the parameter should always be positive
			// make sure that m stays at least 0
			if (*x >= 0 && *m + aux >= 0)
        	{
        		*m = *m >= 0 ? *m : 0;
        	}
        }
        if (*x > *M)
        {
            // make M larger
            *M = *x + (*M - *m) * 0.1;
        }

        printf("----           The range has been adjusted to [%g, %g].\n", *m, *M);

        return (EXIT_SUCCESS);
    }

    // if x is right on the range border, update it to move it away from the border
    if (*x == *m)
    {
        *x = *m + (*M - *m) * 0.000009;
    }
    if (*x == *M)
    {
        *x = *M - (*M - *m) * 0.000009;
    }

    return (EXIT_SUCCESS);
}

void fit_gamma(double len, double *y, double *mean, double *shape)
{
    // y contains observed values from a gamma distribution
    unsigned i = 0;
    (*mean) = 0;
    double var = 0;

    for (i = 0; i < len; i++)
    {
        (*mean) += y[i];
        var += y[i] * y[i];
    }

    (*mean) /= len;
    var = (var - (*mean) * (*mean) * len) / (len - 1.0);
    (*shape) = (*mean) * (*mean) / var;
}

void set_expec_sel_to_obs(Params *p, double **expec)
{
    // calculate the mean observed selected counts
    // and store them in the selected expectations
    unsigned i = 0;
    unsigned j = 0;

    for (j = 0; j < p->pm->n - 1; j++)
    {
        (*expec)[j] = 0;
        for (i = 0; i < p->no_sel; i++)
        {
            (*expec)[j] += p->counts_sel[i].sfs[j] / p->counts_sel[i].len_poly;
        }
        (*expec)[j] /= p->no_sel;
    }

    // divergence counts have a different length - treat them separately
    j = p->pm->n - 1;
    (*expec)[j] = 0;
    for (i = 0; i < p->no_sel; i++)
    {
        (*expec)[j] += p->counts_sel[i].sfs[j] / p->counts_sel[i].len_div;
    }
    (*expec)[j] /= p->no_sel;
}

void set_obs_diff(Params *p, int i, int j, double **obs, double *obs_sel)
{
    int len = p->counts_neut[i].len_poly;
    // divergence counts have a different length
    if (j == p->pm->n)
    {
        len = p->counts_neut[i].len_div;
    }

    // calculate tilde{p_j} and tilde{p_{n-j}}
    (*obs)[0] = p->counts_neut[i].sfs[j - 1] / len;

    if (j < p->pm->n)
    {
        (*obs)[1] = p->counts_neut[i].sfs[p->pm->n - j - 1] / len;
    }
}

int estimate_mut(Params *p)
{
    // from the observed counts of neutral singletons
    // and for given ancestral and contamination errors
    // estimate
    //      the mean and shape of the gamma distribution of the mutation rates
    //      r_{n-1}
    unsigned i = 0;
    double aux = 1 - 2 * p->pm->eps_an;
    double mut_mean = 0;
    double mut_shape = 0;
    double r_n = 0;
    double r_d = 0;
    double *obs_diff = malloc(sizeof(double) * 2);
    // store the different mutation rates estimated from each fragment
    double *obs_theta = malloc(sizeof(double) * p->no_neut);

    for (i = 0; i < p->no_neut; i++)
    {
        set_obs_diff(p, i, 1, &obs_diff, p->expec_sel);
        obs_theta[i] = ((1 - p->pm->eps_an) * obs_diff[0]
                        - p->pm->eps_an * obs_diff[1]) / aux;
        r_n += (1 - p->pm->eps_an) * obs_diff[1] - p->pm->eps_an * obs_diff[0];
        r_d += (1 - p->pm->eps_an) * obs_diff[0] - p->pm->eps_an * obs_diff[1];
    }

    fit_gamma(p->no_neut, obs_theta, &mut_mean, &mut_shape);
    if (p->pm->theta_bar_flag != FALSE)
    {
        p->pm->theta_bar = mut_mean;
        // make sure it's positive
        p->pm->theta_bar = p->pm->theta_bar > 0
        		? p->pm->theta_bar
				: (p->pm->theta_bar_min
						+ (p->pm->theta_bar_max - p->pm->theta_bar_min) * 0.0001);
    }
    if (p->pm->a_flag != FALSE)
    {
    	p->pm->a = mut_shape;
    	// make sure it's positive
    	p->pm->a = p->pm->a > 0
    			? p->pm->a
				: (p->pm->a_min
						+ (p->pm->a_max - p->pm->a_min) * 0.0001);
    }

    free(obs_theta);
    free(obs_diff);

    // for some reason, (*obs_r) can be negative sometimes
    aux = (p->pm->n - 1) * r_n / r_d;
    return (aux > 0 ? aux : 1);
}

void estimate_neut(Params *p)
{
    // from the observed count of neutral polymorphism and divergence
    // and for given ancestral and contamination errors
    // estimate
    //      the mean and shape of the gamma distribution of the mutation rates
    //      r_i (demography)
    //      lambda (divergence)
    double aux = 0;

    if (p->pm->theta_bar_flag != FALSE || p->pm->a_flag != FALSE
                    || p->pm->r_flag != FALSE)
    {
        aux = estimate_mut(p);
    }

    double *obs_diff = malloc(sizeof(double) * 2);
    // store the n-1 observed r parameters for each type of count
    double *obs_r = malloc(sizeof(double) * (p->pm->n - 1));
    unsigned i = 0;
    unsigned j = 0;
    unsigned total = 0;
    unsigned groups = p->pm->n - 1 < p->pm->no_groups ? p->pm->n - 1 : p->pm->no_groups;

    if (p->pm->r_flag != FALSE)
    {
        // calculate all n-1 possible r values
        set_to_zero(&obs_r, p->pm->n - 1);
        // r_{n-1} is calculated from estimate_mut and stored in aux
        obs_r[p->pm->n - 2] = aux;
        aux = p->no_neut * p->pm->theta_bar * (1 - 2 * p->pm->eps_an);

        // skip first and last positions
        // first is always set to 1, last is already calculated
        for (j = 1; j < p->pm->n - 2; j++)
        {
            for (i = 0; i < p->no_neut; i++)
            {
                set_obs_diff(p, i, j + 1, &obs_diff, p->expec_sel);
                obs_r[j] += (1 - p->pm->eps_an) * obs_diff[0]
                                + p->pm->eps_an * obs_diff[1];
            }
            obs_r[j] *= (j + 1) / aux;
        }

        // initialize the rs to 0
        set_to_zero(&p->pm->r, p->pm->no_groups);
        // the first r is always 1
        p->pm->r[0] = 1;

        // set the grouped r by using the means
        j = 1;
        for (i = 1; i < groups; i++)
        {
            aux = 0;
            total = 0;
            while (j < p->pm->n - 1 && p->pm->inv_groups[j] == i)
            {
            	// r values can be negative sometimes
            	// make sure they are 1 in that case
            	obs_r[j] = obs_r[j] > 0 ? obs_r[j] : 1;
                aux += obs_r[j];
                total++;
                j++;
            }
            if (total != 0)
            {
                p->pm->r[i] = aux / total;
            }
            else
            {
                p->pm->r[i] = 0;
            }
        }

        // if divergence r is not grouped, I need to copy r_{n-1} to r_{n}
        if (p->pm->no_r == p->pm->n && p->pm->r[p->pm->no_groups - 1] == 0)
        {
            p->pm->r[p->pm->no_groups - 1] = p->pm->r[p->pm->no_groups - 2];
        }
    }

    if (p->pm->lambda_flag != FALSE)
    {
        for (i = 0; i < p->no_neut; i++)
        {
            set_obs_diff(p, i, p->pm->n, &obs_diff, p->expec_sel);
            p->pm->lambda += obs_diff[0];
        }
        if (p->pm->no_r == p->pm->n)
        {
            // I am using r for divergence
            p->pm->lambda = p->pm->lambda / p->no_neut
                            - p->pm->theta_bar * p->pm->r[p->pm->no_groups - 1]
                                            / p->pm->n;
        }
        else
        {
            p->pm->lambda = p->pm->lambda / p->no_neut
                            - p->pm->theta_bar / p->pm->n;
        }
        // make sure it's positive
        p->pm->lambda = p->pm->lambda > 0
                		? p->pm->lambda
        				: (p->pm->lambda_min
        						+ (p->pm->lambda_max - p->pm->lambda_min) * 0.0001);
    }

    free(obs_diff);
    free(obs_r);
}

int estimate_grid_neut(Params *p, int no_grids)
{
    // do a grid search for eps ancestral and eps contamination
    double best_lk = -DBL_MAX;
    int status = EXIT_SUCCESS;

    ParamsModel best_pm;
	initialize_params_model(&best_pm);
	best_pm.n = p->pm->n;
	best_pm.no_sel = p->pm->no_sel;
	allocate_selection_params(&best_pm);
	allocate_grouping(&best_pm, NULL);
	copy_params_model(&best_pm, *p->pm);

    set_expec_sel_to_obs(p, &p->expec_sel);

    int it = 0;
    int i = 0;

    if (p->pm->eps_an_flag == FALSE)
    {
        estimate_neut(p);
        copy_params_model(&best_pm, *p->pm);
    }
    else
    {
        for (it = 0; it <= no_grids; it++)
        {
            set_params_eps(p->pm, no_grids, it);
            estimate_neut(p);
            set_neut_lnL(p);
            if (p->lnL_neut > best_lk)
            {
                best_lk = p->lnL_neut;
                copy_params_model(&best_pm, *p->pm);
            }
        }
    }

    copy_params_model(p->pm, best_pm);
    free_params_model(&best_pm);

    // check the ranges for the parameters
    if (p->pm->theta_bar_flag != FALSE)
    {
        status += check_lim_update(&p->pm->theta_bar, &p->pm->theta_bar_min,
                                   &p->pm->theta_bar_max, TRUE, "theta_bar");
    }
    if (p->pm->a_flag != FALSE)
    {
        status += check_lim_update(&p->pm->a, &p->pm->a_min, &p->pm->a_max,
                                   TRUE, "a");
    }
    if (p->pm->r_flag != FALSE)
    {
        char *aux = (char*) malloc(20 * sizeof(char));
        for (i = 1; i < p->pm->no_groups; i++)
        {
        	sprintf(aux, "r %d", i + 1);
            status += check_lim_update(&p->pm->r[i], &p->pm->r_min,
                                       &p->pm->r_max, TRUE, aux);
        }
        free(aux);
    }
    if (p->pm->lambda_flag != FALSE)
    {
        status += check_lim_update(&p->pm->lambda, &p->pm->lambda_min,
                                   &p->pm->lambda_max, TRUE, "lambda");
    }

    return (status);
}

void estimate_grid_sel(Params *p, int no_grids, int no_params)
{
    if (p->pm->model != 4)
    {
        // do a grid search for selection parameters
        double best_lk = -DBL_MAX;

        ParamsModel best_pm;
        initialize_params_model(&best_pm);
        best_pm.n = p->pm->n;
        best_pm.no_sel = p->pm->no_sel;
        allocate_selection_params(&best_pm);
        allocate_grouping(&best_pm, NULL);
        copy_params_model(&best_pm, *p->pm);

        int it[4] = { 0, 0, 0, 0 };
        for (it[0] = 0; it[0] <= no_grids; it[0]++)
        {
            for (it[1] = 0; it[1] <= no_grids; it[1]++)
            {
                for (it[2] = 0; it[2] <= no_grids; it[2]++)
                {
                    for (it[3] = 0; it[3] <= no_grids; it[3]++)
                    {
                        // set parameters
                        set_params_sel(p->pm, no_grids, it);
                        set_sel_lnL(p);
                        // s_bar has to be smaller than s_max
                        // otherwise, I get Inf
                        if (p->lnL_sel != DBL_MAX && p->lnL_sel > best_lk)
                        {
                            best_lk = p->lnL_sel;
                            copy_params_model(&best_pm, *p->pm);
                        }
                        if (no_params < 4)
                        {
                            break;
                        }
                    }
                    if (no_params < 3)
                    {
                        break;
                    }
                }
                if (no_params < 2)
                {
                    break;
                }
            }
        }
        copy_params_model(p->pm, best_pm);
        free_params_model(&best_pm);
    }
    else
    {
        // set everything to uniform for now
        // TODO I will need to update this to something better
        unsigned i = 0;
        double prob = 1.0 / (p->pm->no_sel / 2);
        for (i = 0; i < p->pm->no_sel/2; i++)
        {
            p->pm->sel_params[2*i+1] = prob;
        }
    }
}

int estimate_grid(ParamsShare *ps)
{
	int status = EXIT_SUCCESS;
    int no_grids = 10;
    size_t i = 0;

    for (i = 0; i < ps->no_data; i++)
    {
        if (ps->p[i].pm->neut > 0 || ps->neut > 0)
        {
            status += estimate_grid_neut(&ps->p[i], no_grids);
        }

        if (ps->p[i].pm->sel > 0 || ps->sel > 0)
        {
            // if I am not using model A, set a lower number of grids
            if (ps->p[i].pm->model != 1)
            {
                no_grids = 6;
            }
            // have to specify how many parameters I am estimating
            estimate_grid_sel(&ps->p[i], no_grids, ps->sel + ps->p[i].pm->sel);
        }
    }

    if (status != EXIT_SUCCESS)
    {
        fprintf(stderr,
                "---- Warning: some of the estimated parameters did not fit "
                "the given ranges. The ranges have been adjusted automatically.\n");
    }

    return (EXIT_SUCCESS);
}

void initialize_shared(ParamsShare *ps)
{
    unsigned i = 0;
    unsigned j = 0;
    double value = 0;
    double range_min = 0;
    double range_max = 0;

    // the initial value of the shared parameters
    // is given by the average value of the
    // of the initial value from Params

    // the range of the shared parameters
    // is also given by the average ranges from Params

    if (ps->use_neut_ln == TRUE)
    {
        if (ps->r_shared == TRUE)
        {
            // r values for the average
            double *r = malloc(sizeof(double) * ps->no_groups);
            range_min = 0;
            range_max = 0;
            // counter for how many files have each r value
            double *count = malloc(sizeof(double) * ps->no_groups);
            set_to_zero(&r, ps->no_groups);
            set_to_zero(&count, ps->no_groups);
            // get the sum
            for (j = 0; j < ps->no_data; j++)
            {
                for (i = 1; i < ps->p[j].pm->no_groups; i++)
                {
                    r[i] += ps->p[j].pm->r[i];
                    count[i]++;
                }
                range_min += ps->p[j].pm->r_min;
                range_max += ps->p[j].pm->r_max;
            }
            // get the average
            for (i = 1; i < ps->no_groups; i++)
            {
                r[i] /= count[i];
            }
            range_min /= ps->no_data;
            range_max /= ps->no_data;
            // copy back the values
            for (j = 0; j < ps->no_data; j++)
            {
                for (i = 1; i < ps->p[j].pm->no_groups; i++)
                {
                    ps->p[j].pm->r[i] = r[i];
                }
                ps->p[j].pm->r_min = range_min;
                ps->p[j].pm->r_max = range_max;
            }
            // free memory
            free(r);
            free(count);
        }
        if (ps->lambda_shared == TRUE)
        {
            value = 0;
            range_min = 0;
            range_max = 0;
            // calculate mean
            for (j = 0; j < ps->no_data; j++)
            {
                value += ps->p[j].pm->lambda;
                range_min += ps->p[j].pm->lambda_min;
                range_max += ps->p[j].pm->lambda_max;
            }
            value /= ps->no_data;
            range_min /= ps->no_data;
            range_max /= ps->no_data;
            // copy back
            for (j = 0; j < ps->no_data; j++)
            {
                ps->p[j].pm->lambda = value;
                ps->p[j].pm->lambda_min = range_min;
                ps->p[j].pm->lambda_max = range_max;
            }
        }
        if (ps->theta_bar_shared == TRUE)
        {
            value = 0;
            range_min = 0;
            range_max = 0;
            // calculate mean
            for (j = 0; j < ps->no_data; j++)
            {
                value += ps->p[j].pm->theta_bar;
                range_min += ps->p[j].pm->theta_bar_min;
                range_max += ps->p[j].pm->theta_bar_max;
            }
            value /= ps->no_data;
            range_min /= ps->no_data;
            range_max /= ps->no_data;
            // copy back
            for (j = 0; j < ps->no_data; j++)
            {
                ps->p[j].pm->theta_bar = value;
                ps->p[j].pm->theta_bar_min = range_min;
                ps->p[j].pm->theta_bar_max = range_max;
            }
        }
        if (ps->a_shared == TRUE)
        {
            value = 0;
            range_min = 0;
            range_max = 0;
            // calculate mean
            for (j = 0; j < ps->no_data; j++)
            {
                value += ps->p[j].pm->a;
                range_min += ps->p[j].pm->a_min;
                range_max += ps->p[j].pm->a_max;
            }
            value /= ps->no_data;
            range_min /= ps->no_data;
            range_max /= ps->no_data;
            // copy back
            for (j = 0; j < ps->no_data; j++)
            {
                ps->p[j].pm->a = value;
                ps->p[j].pm->a_min = range_min;
                ps->p[j].pm->a_max = range_max;
            }
        }
        if (ps->eps_an_shared == TRUE)
        {
            value = 0;
            range_min = 0;
            range_max = 0;
            // calculate mean
            for (j = 0; j < ps->no_data; j++)
            {
                value += ps->p[j].pm->eps_an;
                range_min += ps->p[j].pm->eps_an_min;
                range_max += ps->p[j].pm->eps_an_max;
            }
            value /= ps->no_data;
            range_min /= ps->no_data;
            range_max /= ps->no_data;
            // copy back
            for (j = 0; j < ps->no_data; j++)
            {
                ps->p[j].pm->eps_an = value;
                ps->p[j].pm->eps_an_min = range_min;
                ps->p[j].pm->eps_an_max = range_max;
            }
        }
    }

    if (ps->use_sel_ln == TRUE)
    {
        // I don't need to treat model D separately because I have the flags set
        // to FALSE for the selection coefficients
        for (i = 0; i < ps->p[0].pm->no_sel; i++)
        {
            if (ps->sel_shared[i] == TRUE)
            {
                value = 0;
                range_min = 0;
                range_max = 0;
                // calculate mean
                for (j = 0; j < ps->no_data; j++)
                {
                    value += ps->p[j].pm->sel_params[i];
                    range_min += ps->p[j].pm->sel_min[i];
                    range_max += ps->p[j].pm->sel_max[i];
                }
                value /= ps->no_data;
                range_min /= ps->no_data;
                range_max /= ps->no_data;
                // copy back
                for (j = 0; j < ps->no_data; j++)
                {
                    ps->p[j].pm->sel_params[i] = value;
                    ps->p[j].pm->sel_min[i] = range_min;
                    ps->p[j].pm->sel_max[i] = range_max;
                }
            }
        }
    }
}

/****************************************************************************
 * Functions on structures:
 * initialization, allocation, freeing and printing
 ****************************************************************************/
void initialize_params_basin_hop(ParamsBasinHop *pb)
{
    pb->max_iter = 0;

    pb->max_same = 50;
    pb->temp = 1;
    pb->step = 50;
    pb->accept_rate = 0.5;
    pb->interval = 10;
    pb->factor = 0.9;

    pb->no_step = 0;
    pb->no_accept = 0;

    pb->size = 0;
    pb->x = NULL;

    pb->lnL = DBL_MAX;
    pb->criteria = DBL_MAX;
    pb->best_lnL = DBL_MAX;

    // initialize random number generator
    srand(time(NULL));
}

void set_default_max_iter(ParamsBasinHop *pb)
{
    pb->max_iter = 500;
}

void allocate_params_basin_hop(ParamsBasinHop *pb, ParamsShare ps)
{
    pb->size = ps.neut + ps.sel;

    if (pb->size > 0)
    {
        pb->x = gsl_vector_alloc(pb->size);
    }

    // setup random number generator
    srand(time(NULL));
}

void free_params_basin_hop(ParamsBasinHop *pb)
{
    if (pb->x)
    {
        gsl_vector_free(pb->x);
    }
}

/****************************************************************************
 * Functions used for storing best optimum
 ****************************************************************************/
int add_best_optimum(ParamsBasinHop *pb, ParamsShare *ps, int gradient_descent)
{
    double thresh = 0.25;

    // I compare the lk of the new optimum to the all time best
    // I accept new optimum if
    // the new solution is a lot better, p->lnL < best_lnL - thresh
    // or the new solution is better or just slightly worse
    // but the gradient of the new solution is better
    // p->lnL < best_lnL + thresh && p->grad < pb->grad
    if ((ps->lnL < pb->best_lnL - thresh)
        | ((ps->lnL < pb->best_lnL + thresh) & (ps->criteria < pb->criteria)))
    {
        // update the best likelihood ever found
        pb->best_lnL = ps->lnL < pb->best_lnL ? ps->lnL : pb->best_lnL;

        // update best solution
        pb->lnL = ps->lnL;
        pb->criteria = ps->criteria;
        transform(&pb->x, *ps);

        if (ps->verbose > 0)
        {
            printf("-- Found new optimum with likelihood %.15f ", -ps->lnL);
            if (gradient_descent == TRUE)
            {
            	printf("and gradient %.5f\n\n", ps->criteria);
            }
            else
            {
            	printf("and size %.5f\n\n", ps->criteria);
            }
        }

        return (TRUE);
    }
    printf("\n");

    // update the best likelihood ever found
    pb->best_lnL = ps->lnL < pb->best_lnL ? ps->lnL : pb->best_lnL;

    return (FALSE);
}

/****************************************************************************
 * Metropolis Hastings acceptance criterion
 ****************************************************************************/
int accept_reject(ParamsBasinHop *pb, double lnL_new, double lnL_old)
{
    double w = exp(-(lnL_new - lnL_old) / pb->temp);
    double u =  rand() / RAND_MAX;

    if (w >= u)
    {
        pb->no_accept++;
        return (TRUE);
    }

    return (FALSE);
}

/****************************************************************************
 * Adaptive step taking
 ****************************************************************************/
void take_adaptive_step(ParamsBasinHop *pb, ParamsShare *ps, gsl_vector **x)
{
    undo_transform(ps, pb->x);
    random_step(ps, pb->step);
    transform(x, *ps);

    pb->no_step++;
    if (pb->no_step % pb->interval == 0)
    {
        // adjust step size
        double accept_rate = ((double) pb->no_accept) / pb->no_step;
        if (accept_rate > pb->accept_rate)
        {
            // I accept too many steps. Take bigger steps
            // this generally means I am trapped in a basin
            pb->step /= pb->factor;
            if (ps->verbose > 0)
            {
                printf("---- Current accepting rate %g is too large, "
                       "increasing step size to %g\n",
                       accept_rate, pb->step);
            }
        }
        else
        {
            // I accept too few steps. Take smaller steps.
            pb->step *= pb->factor;
            if (ps->verbose > 0)
            {
                printf("---- Current accepting rate %g is too small, "
                       "decreasing step size to %g\n",
                       accept_rate, pb->step);
            }
        }
    }
}

void print_final_result(ParamsShare *ps)
{
    unsigned i = 0;
    unsigned j = 0;

    double curr_lambda = 0;
    double mis_neut = 0;
    double mis_sel_full = 0;
    double mis_sel_del = 0;
    double expec_fixed_del = 0;
    double expec_fixed_full = 0;
    double div_neut = 0;
    double div_sel = 0;
    double len_neut = 0;
    double len_sel = 0;
    unsigned div_flag = 0;

    for (j = 0; j < ps->no_data; j++)
    {
        printf("\n---- Results for %s \n", ps->data_files[j]);

        fprintf_params_model(*ps->p[j].pm, stdout, "-- ");
        printf("\n");

		// reset numbers to 0
		mis_neut = 0;
		mis_sel_full = 0;
		mis_sel_del = 0;
		expec_fixed_del = 0;
		expec_fixed_full = 0;
		div_neut = 0;
		div_sel = 0;
		len_neut = 0;
		len_sel = 0;

        // print expected SFS, both neutral and selected, per site
        printf("---- Expected P_neut(i), 0 < i < n (neutral SFS per site) \n");
        for (i = 0; i < ps->p[j].pm->n - 1; i++)
        {
            printf("E[P_neut(%d)] = %10.10f\n", i + 1, ps->p[j].expec_neut[i]);
        }
        printf("\n---- Expected P_sel(i), 0 < i < n (selected SFS per site) \n");
        for (i = 0; i < ps->p[j].pm->n - 1; i++)
        {
            printf("E[P_sel(%d)] = %10.10f\n", i + 1, ps->p[j].expec_sel[i]);
        }

        // calculate observed divergence counts
        for (i = 0; i < ps->p[j].no_neut; i++)
        {
            div_neut += ps->p[j].counts_neut[i].sfs[ps->p[j].pm->n - 1];
            len_neut += ps->p[j].counts_neut[i].len_div;
        }
        for (i = 0; i < ps->p[j].no_sel; i++)
        {
            div_sel += ps->p[j].counts_sel[i].sfs[ps->p[j].pm->n - 1];
            len_sel += ps->p[j].counts_sel[i].len_div;
        }

        if (ps->p[j].pm->div_flag == TRUE)
        {
            // print expected divergence counts
            printf("\n---- Expected D_neut and D_sel "
                   "(neutral and selected divergence per site) \n");
            printf("E[D_neut] = %10.10f\nE[D_sel] = %10.10f\n",
                   ps->p[j].expec_neut[ps->p[j].pm->n - 1],
                   ps->p[j].expec_sel[ps->p[j].pm->n - 1]);
        }

        div_flag = ps->p[j].pm->div_flag;

        // for the rest, I need p->pm->div_flag = TRUE
        ps->p[j].pm->div_flag = TRUE;

        // in preparation for alpha, calculate the integral over the del DFE from eq. 8
        // and the integral over the full DFE from eq. 19
        // I can obtain this by getting the expected divergence with lambda = 1
        // and then later on subtract the misattributed polymorphism
        curr_lambda = ps->p[j].pm->lambda;
        ps->p[j].pm->lambda = 1;
        set_sel_expec(ps->p[j].pm, &ps->p[j].expec_sel, FALSE);
        expec_fixed_full = ps->p[j].expec_sel[ps->p[j].pm->n - 1];
        set_sel_expec(ps->p[j].pm, &ps->p[j].expec_sel, TRUE);
        expec_fixed_del = ps->p[j].expec_sel[ps->p[j].pm->n - 1];
        ps->p[j].pm->lambda = curr_lambda;

        // calculate missatribuated polymorphism
        // I can do this easily by setting lambda to 0
        curr_lambda = ps->p[j].pm->lambda;
        ps->p[j].pm->lambda = 0;
        set_sel_expec(ps->p[j].pm, &ps->p[j].expec_sel, FALSE);
        mis_sel_full = ps->p[j].expec_sel[ps->p[j].pm->n - 1];
        set_sel_expec(ps->p[j].pm, &ps->p[j].expec_sel, TRUE);
        mis_sel_del = ps->p[j].expec_sel[ps->p[j].pm->n - 1];
        set_neut_expec(*ps->p[j].pm, &ps->p[j].expec_neut);
        mis_neut = ps->p[j].expec_neut[ps->p[j].pm->n - 1];

        // if I am not using divergence data, the below messages are pointless
        if (ps->p[j].counts_neut[0].sfs[ps->p[j].pm->n - 1] != -1)
        {
            printf("\n---- Expected neutral and selected misattributed "
                   "polymorphism per site\n");
            if (div_flag == FALSE)
            {
                fprintf(stderr,
                        "Option -w was used: r_n was not estimated, and it is set to 1 "
                        "(i.e. no demography effect) while calculating misattributed polymorphism.\n");
            }
            printf("E[mis_neut] = %10.10f\nE[mis_sel] = %10.10f\n", mis_neut,
                   mis_sel_full);
        }

        // reset lambda and re-calculate expectations
        ps->p[j].pm->lambda = curr_lambda;
        set_sel_lnL_share(ps);
        set_neut_lnL_share(ps);

        // remove misattributed polymorphism from expectations of fixed mutations
        expec_fixed_full -= mis_sel_full;
        expec_fixed_del -= mis_sel_del;

        // if I am not using divergence data, alpha_div should not be printed
        if (ps->p[j].counts_neut[0].sfs[ps->p[j].pm->n - 1] != -1)
        {
            // re-adjust divergence counts for the misattributed polymorhphism (eq. 9)
            div_neut -= len_neut * mis_neut;
            div_sel -= len_sel * mis_sel_full;
            // calculate alpha_div (eq. 8)
            printf("\n---- alpha_div = %.6f",
                   1 - len_sel * div_neut * expec_fixed_del / (len_neut * div_sel));
        }

        // calculate alpha_dfe (eq. 10)
        printf("\n---- alpha_dfe = %.6f",
                       1 - expec_fixed_del / expec_fixed_full);

        printf("\n\n");
    }
}

/****************************************************************************
 * The basin hopping algorithm
 ****************************************************************************/
int run_basin_hopping(ParamsBasinHop *pb,
                      const gsl_multimin_fdfminimizer_type *type,
                      ParamsOptim po, ParamsShare *ps)
{
    int status = EXIT_SUCCESS;

    gsl_vector *x = gsl_vector_alloc(pb->size);

    // char array for printing purposes
    char *s = malloc(sizeof(char) * 50);

    unsigned iter = 0;
    unsigned same_cnt = 0;

    double lnL = DBL_MAX;
    int new_min = FALSE;

    double max_iter = 0;
    int init = 0;

    // do I want to run both initial estimation and automatic estimation?
    // when min_init == 0 and max_init == 1, it runs both
    int min_init = 0;
    int max_init = 0;
    if (ps->initial_values == FALSE && ps->initial_estimation == TRUE)
    {
        // only automatic estimation
        min_init = 1;
        max_init = 1;
    }
    if (ps->initial_values == TRUE && ps->initial_estimation == TRUE)
    {
        // both
        max_init = 1;
    }

    // flag for final printing - it is only set to TRUE if it doesn't fail with errors
    unsigned print_flag = FALSE;

    // make sure it works both on neutral and selected
    ps->use_neut_ln = TRUE;
    ps->use_sel_ln = TRUE;

    if (pb->size > 0)
    {
        if (pb->max_iter > 0 && pb->size > 0)
        {
            printf("\n---- Starting a maximum of %.0f basin hopping iterations\n",
                   pb->max_iter);
        }
        for (init = min_init; init <= max_init; init++)
        {
            ps->initial_estimation = init;

            // initialize parameters
            if (ps->initial_estimation == TRUE)
            {
                printf("\n---- Calculating initial values\n");
                status = estimate_grid(ps);
            }
            else
            {
                printf("\n---- Using provided initial values\n");
            }

            // have to make sure that the shared parameters
            // are initialized as the mean over Params
            initialize_shared(ps);

            // initialize local minimum
            add_best_optimum(pb, ps, po.grad_descent);

            // initialize likelihood
            set_sel_lnL_share(ps);
            set_neut_lnL_share(ps);
            lnL = ps->lnL_neut + ps->lnL_sel;
            if (lnL == DBL_MAX)
            {
                status = FOUND_NAN;
            }

            // initialize likelihood to DBL_MAX
            // so that the rest of basin hopping works
            lnL = DBL_MAX;

            // initialize vector for storing intermediate solutions
            transform(&x, *ps);

            if (status != FOUND_NAN)
            {
                // when running basin hopping and min_init != max_init
                // I want to run basin hopping on max_init only
                max_iter = init == max_init ? pb->max_iter : 0;
                // run the basin hopping iterations
                for (iter = 0; iter <= max_iter; iter++)
                {
                    printf("-- Starting local optimization\n");

                    status = optimize(type, po, ps, s);
                    if (status != EXIT_SUCCESS)
                    {
                        break;
                    }

                    // accept the move based on Metropolis Hastings test
                    int accept = accept_reject(pb, ps->lnL, lnL);
                    if (accept == TRUE)
                    {
                        if (ps->verbose > 0 && max_iter > 0)
                        {
                            printf("-- Accepted new set of parameters with "
                                   "likelihood %.15f", -ps->lnL);
                            if (po.grad_descent == TRUE)
                            {
                            	printf(" and gradient %.5f\n", ps->criteria);
                            }
                            else
                            {
                            	printf(" and size %.5f\n", ps->criteria);
                            }
                        }
                        new_min = add_best_optimum(pb, ps, po.grad_descent);
                        // update current solution
                        lnL = ps->lnL;
                        transform(&x, *ps);
                    }
                    else
                    {
                        printf("\n");
                    }

                    // test if I am stuck in the same solution
                    same_cnt++;
                    if (new_min == TRUE)
                    {
                        same_cnt = 0;
                    }
                    else if (same_cnt >= pb->max_same)
                    {
                        status = SAME_STATE;
                        break;
                    }

                    if (1 <= max_iter && iter < max_iter)
                    {
                        printf("---- Basin hopping iteration %d\n", iter + 1);
                    }

                    // randomly displace the coordinates
                    take_adaptive_step(pb, ps, &x);
                }

                // test if reached max iterations
                status = iter >= max_iter ? MAX_ITER : status;
                if (max_iter > 0 && pb->size > 0)
                {
                    printf("\n---- Basin hopping performed %d iterations\n",
                           iter - 1);
                    switch (status)
                    {
                        case MAX_ITER:
                        {
                            printf("---- Basin hopping reached maximum number of "
                                            "iterations allowed\n\n");
                            break;
                        }
                        case SAME_STATE:
                        {
                            printf("---- Basin hopping is stuck in the same solution\n\n");
                            break;
                        }
                        case GSL_SUCCESS:
                        {
                            printf("---- Basin hopping found a good candidate solution\n\n");
                            break;
                        }
                    }
                }

                // make sure the parameters contain the best found solution
                undo_transform(ps, pb->x);
            }
            else
            {
                fprintf(stderr, "\n---- Starting parameters have nan likelihood\n");
                fprintf(stderr, "\n");
            }
        }

        // print best solution found
        if (status != FOUND_NAN && pb->lnL != DBL_MAX)
        {

            printf("---- Best joint likelihood found %.15f", -pb->lnL);
            if (po.grad_descent == TRUE)
            {
            	printf(" with gradient %.5f\n", pb->criteria);
            }
            else
            {
            	printf(" with size %.5f\n", pb->criteria);
            }
            // set final print flag
            print_flag = TRUE;
            // make sure I have the correct expectations in p
            set_sel_lnL_share(ps);
            set_neut_lnL_share(ps);
        }
    }
    else
    {
        // nothing to optimize, calculate joint likelihood
        set_sel_lnL_share(ps);
        set_neut_lnL_share(ps);
        lnL = -(ps->lnL_neut + ps->lnL_sel);
        pb->lnL = lnL;
        printf("\n---- Joint likelihood %.15f\n", -lnL);
        // set final print flag
        print_flag = TRUE;
    }

    if (print_flag == TRUE && status != FOUND_NAN && pb->lnL != DBL_MAX)
    {
        print_final_result(ps);
    }

    // free memory
    gsl_vector_free(x);
    free(s);

    return (EXIT_SUCCESS);
}
