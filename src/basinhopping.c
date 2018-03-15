/*
 * polyDFE v1.0: predicting DFE and alpha from polymorphism data
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
    if (flag == FALSE)
    {
        return (EXIT_SUCCESS);
    }

    // verify that m <= x <= M
    if (*x < *m || *x > *M)
    {
        // update limits so that x is within the limits
        if (*x < *m)
        {
            // make m smaller
            // this is only checked for parameters that are always positive
            // so make sure that m stays at least 0
            *m = *x - (*M - *m) * 0.1;
            *m = *m >= 0 ? *m : 0;
        }
        if (*x > *M)
        {
            // make M larger
            *M = *x + (*M - *m) * 0.1;
        }

        printf("---- The initial value of %s did not fit the given ranges. "
               "The ranges have been adjusted automatically.\n", s);

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

int estimate_mut(Params *p, double *obs_r)
{
    // from the observed counts of neutral singletons
    // and for given ancestral and contamination errors
    // estimate
    //      the mean and shape of the gamma distribution of the mutation rates
    //      r_{n-1}
    int status = EXIT_SUCCESS;
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
    if (p->pm->theta_bar_flag == TRUE)
    {
        p->pm->theta_bar = mut_mean;
        // check estimation is within given ranges
        status += check_lim_update(&p->pm->theta_bar, &p->pm->theta_bar_min,
                                   &p->pm->theta_bar_max, TRUE, "theta_bar");
    }
    if (p->pm->a_flag == TRUE)
    {
        p->pm->a = mut_shape;
        // check estimation is within given ranges
        status += check_lim_update(&p->pm->a, &p->pm->a_min, &p->pm->a_max,
                                   TRUE, "a");
    }

    // for some reason, (*obs_r) can be negative sometimes
    aux = (p->pm->n - 1) * r_n / r_d;
    (*obs_r) = aux > 0 ? aux : 1;

    free(obs_theta);
    free(obs_diff);

    return (status);
}

int estimate_neut(Params *p)
{
    // from the observed count of neutral polymorphism and divergence
    // and for given ancestral and contamination errors
    // estimate
    //      the mean and shape of the gamma distribution of the mutation rates
    //      r_i (demography)
    //      lambda (divergence)
    int status = EXIT_SUCCESS;
    double aux = 0;

    if (p->pm->theta_bar_flag == TRUE || p->pm->a_flag == TRUE
                    || p->pm->r_flag == TRUE)
    {
        status = estimate_mut(p, &aux);
        if (status != EXIT_SUCCESS)
        {
            return (status);
        }
    }

    double *obs_diff = malloc(sizeof(double) * 2);
    // store the n-1 observed r parameters for each type of count
    double *obs_r = malloc(sizeof(double) * (p->pm->n - 1));
    unsigned i = 0;
    unsigned j = 0;
    unsigned total = 0;
    unsigned groups = p->pm->n - 1 < p->pm->no_groups ? p->pm->n - 1 : p->pm->no_groups;

    if (p->pm->r_flag == TRUE)
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

        // check estimation is within given ranges
        for (i = 1; i < p->pm->no_groups; i++)
        {
            status += check_lim_update(&p->pm->r[i], &p->pm->r_min,
                                       &p->pm->r_max, TRUE, "r");
        }
    }

    if (p->pm->lambda_flag == TRUE)
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

        // check estimation is within given ranges
        status += check_lim_update(&p->pm->lambda, &p->pm->lambda_min,
                                   &p->pm->lambda_max, TRUE, "lambda");
    }

    free(obs_diff);
    free(obs_r);

    return (status);
}

int estimate_grid_neut(Params *p, int no_grids)
{
    // do a grid search for eps ancestral and eps contamination
    double best_lk = -DBL_MAX;
    int best_status = EXIT_SUCCESS;
    int curr_status = EXIT_SUCCESS;

    ParamsModel best_pm;
    initialize_params_model(&best_pm);
    best_pm.n = p->pm->n;
    best_pm.no_sel = p->pm->no_sel;
    allocate_selection_params(&best_pm);
    allocate_grouping(&best_pm, NULL);

    set_expec_sel_to_obs(p, &p->expec_sel);

    int it = 0;

    if (p->pm->eps_an_flag == FALSE)
    {
        best_status = estimate_neut(p);
        copy_params_model(&best_pm, *p->pm);
    }
    else
    {
        for (it = 0; it <= no_grids; it++)
        {
            curr_status = EXIT_SUCCESS;
            set_params_eps(p->pm, no_grids, it);
            curr_status += estimate_neut(p);
            set_neut_lnL(p);
            if (p->lnL_neut > best_lk)
            {
                best_lk = p->lnL_neut;
                best_status = curr_status;
                copy_params_model(&best_pm, *p->pm);
            }
        }
    }

    copy_params_model(p->pm, best_pm);
    free_params_model(&best_pm);

    return (best_status);
}

void estimate_grid_sel(Params *p, int no_grids)
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
                        if (p->pm->sel < 4)
                        {
                            break;
                        }
                    }
                    if (p->pm->sel < 3)
                    {
                        break;
                    }
                }
                if (p->pm->sel < 2)
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

int estimate_grid(Params *p)
{
    int status = EXIT_SUCCESS;
    int no_grids = 10;

    if (p->pm->neut > 0)
    {
        printf("---- Calculating initial values for neutral parameters\n");
        status = estimate_grid_neut(p, no_grids);

        if (status != EXIT_SUCCESS)
        {
            fprintf(stderr,
                    "---- Warning: some of the estimated parameters did not fit "
                    "the given ranges. The ranges have been adjusted automatically.\n");
        }
    }

    if (p->pm->sel > 0)
    {
        printf("---- Calculating initial values for selection parameters\n");
        // if I am not using model A, set a lower number of grids
        if (p->pm->model != 1)
        {
            no_grids = 6;
        }
        estimate_grid_sel(p, no_grids);
    }

    return (EXIT_SUCCESS);
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
    pb->grad = DBL_MAX;
    pb->best_lnL = DBL_MAX;

    // initialize random number generator
    srand48(time(NULL));
}

void set_default_max_iter(ParamsBasinHop *pb)
{
    pb->max_iter = 500;
}

void allocate_params_basin_hop(ParamsBasinHop *pb, ParamsModel* pm)
{
    count_params(pm);
    pb->size = pm->neut + pm->sel;

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
int add_best_optimum(ParamsBasinHop *pb, Params *p)
{
    double thresh = 0.25;

    // I compare the lk of the new optimum to the all time best
    // I accept new optimum if
    // the new solution is a lot better, p->lnL < best_lnL - thresh
    // or the new solution is better or just slightly worse
    // but the gradient of the new solution is better
    // p->lnL < best_lnL + thresh && p->grad < pb->grad
    if ((p->lnL < pb->best_lnL - thresh)
        | ((p->lnL < pb->best_lnL + thresh) & (p->grad < pb->grad)))
    {
        // update the best likelihood ever found
        pb->best_lnL = p->lnL < pb->best_lnL ? p->lnL : pb->best_lnL;

        // update best solution
        pb->lnL = p->lnL;
        pb->grad = p->grad;
        transform(&pb->x, *p->pm);

        if (p->verbose > 0)
        {
            printf("-- Found new optimum with likelihood %.15f and gradient "
                            "%.5f\n", -p->lnL, p->grad);
        }
        printf("\n\n");

        return (TRUE);
    }
    printf("\n\n");

    // update the best likelihood ever found
    pb->best_lnL = p->lnL < pb->best_lnL ? p->lnL : pb->best_lnL;

    return (FALSE);
}

/****************************************************************************
 * Metropolis Hastings acceptance criterion
 ****************************************************************************/
int accept_reject(ParamsBasinHop *pb, double lnL_new, double lnL_old)
{
    double w = exp(-(lnL_new - lnL_old) / pb->temp);
    double u =  drand48();

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
void take_adaptive_step(ParamsBasinHop *pb, Params *p, gsl_vector **x)
{
    undo_transform(p->pm, pb->x);
    random_step(p->pm, pb->step);
    transform(x, *p->pm);

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
            if (p->verbose > 0)
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
            if (p->verbose > 0)
            {
                printf("---- Current accepting rate %g is too small, "
                       "decreasing step size to %g\n",
                       accept_rate, pb->step);
            }
        }
    }
}

/****************************************************************************
 * The basin hopping algorithm
 ****************************************************************************/
int run_basin_hopping(ParamsBasinHop *pb,
                      const gsl_multimin_fdfminimizer_type *type,
                      ParamsOptim po, Params *p)
{
    int status = EXIT_SUCCESS;

    gsl_vector *x = gsl_vector_alloc(pb->size);

    // char array for printing purposes
    char *s = malloc(sizeof(char) * 15);

    unsigned iter = 0;
    unsigned same_cnt = 0;

    // these are needed for the prints
    int i = 0;
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

    double lnL = DBL_MAX;
    int new_min = FALSE;

    double max_iter = 0;
    int init = 0;

    // do I want to run both initial estimation and automatic estimation?
    // when min_init == 0 and max_init == 1, it runs both
    int min_init = 0;
    int max_init = 0;
    if (p->pm->initial_values == FALSE && p->pm->inital_estimation == TRUE)
    {
        // only automatic estimation
        min_init = 1;
        max_init = 1;
    }
    if (p->pm->initial_values == TRUE && p->pm->inital_estimation == TRUE)
    {
        // both
        max_init = 1;
    }

    // flag for final printing - it is only set to TRUE if it doesn't fail with errors
    unsigned print_flag = FALSE;

    // check that I have at least two fragments for estimating a
    if (p->no_neut == 1 && p->no_sel == 1 && p->pm->a != -1)
    {
        fprintf(stderr,
                "---- Warning: mutation variability is not used "
                "when only one neutral and one selected fragment is available.\n");
        p->pm->a = -1;
    }

    if (p->pm->a == -1)
    {
        printf("---- No mutation variability. Using Poisson likelihood.\n\n");
        p->pm->a_flag = FALSE;
        p->pm->a_max = 1;
        p->pm->a_min = -2;
    }

    // make sure it works both on neutral and selected
    p->pm->neut_ln = TRUE;
    p->pm->sel_ln = TRUE;

    if (pb->size > 0)
    {
        if (pb->max_iter > 0 && pb->size > 0)
        {
            printf("---- Starting a maximum of %.0f basin hopping iterations\n\n",
                   pb->max_iter);
        }
        for (init = min_init; init <= max_init; init++)
        {
            p->pm->inital_estimation = init;

            // initialize parameters
            if (p->pm->inital_estimation == TRUE)
            {
                status = estimate_grid(p);
            }
            else
            {
                printf("---- Using provided initial values\n");
            }

            // initialize local minimum
            add_best_optimum(pb, p);

            // initialize likelihood
            set_sel_lnL(p);
            set_neut_lnL(p);
            lnL = -(p->lnL_neut + p->lnL_sel);
            if (lnL == - DBL_MAX)
            {
                status = FOUND_NAN;
            }

            // initialize vector for storing intermediate solutions
            transform(&x, *p->pm);

            if (status != FOUND_NAN)
            {
                // when running basin hopping and min_init != max_init
                // I want to run basin hopping on max_init only
                max_iter = init == max_init ? pb->max_iter : 0;
                // run the basin hopping iterations
                for (iter = 0; iter <= max_iter; iter++)
                {
                    printf("-- Starting local optimization\n");

                    status = optimize_using_derivatives(type, po, p, s);
                    if (status != EXIT_SUCCESS)
                    {
                        break;
                    }

                    // accept the move based on Metropolis Hastings test
                    int accept = accept_reject(pb, p->lnL, lnL);
                    if (accept == TRUE)
                    {
                        if (p->verbose > 0 && max_iter > 0)
                        {
                            printf("-- Accepted new set of parameters with "
                                   "likelihood %.15f and gradient "
                                   "%.5f\n",
                                   -p->lnL, p->grad);
                        }
                        new_min = add_best_optimum(pb, p);
                        // update current solution
                        lnL = p->lnL;
                        transform(&x, *p->pm);
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
                    take_adaptive_step(pb, p, &x);
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
                undo_transform(p->pm, pb->x);
            }
            else
            {
                fprintf(stderr, "\n---- Starting parameters have nan likelihood\n");
                fprintf_params_model(*p->pm, stderr, "---- ");
                fprintf(stderr, "\n\n");
            }
        }

        // print best solution found
        if (status != FOUND_NAN && pb->lnL != DBL_MAX)
        {

            printf("---- Best joint likelihood found %.15f with gradient %.5f\n",
                   -pb->lnL, pb->grad);
            fprintf_params_model(*p->pm, stdout, "-- ");
            printf("\n");
            // set final print flag
            print_flag = TRUE;
            // make sure I have the correct expectations in p
            set_sel_lnL(p);
            set_neut_lnL(p);
        }
    }
    else
    {
        // nothing to optimize, calculate joint likelihood
        set_sel_lnL(p);
        set_neut_lnL(p);
        lnL = -(p->lnL_neut + p->lnL_sel);
        pb->lnL = lnL;
        printf("\n---- Joint likelihood %.15f\n", -lnL);
        fprintf_params_model(*p->pm, stdout, "-- ");
        printf("\n");
        // set final print flag
        print_flag = TRUE;
    }

    if (print_flag == TRUE && status != FOUND_NAN && pb->lnL != DBL_MAX)
    {
        // print expected SFS, both neutral and selected, per site
        printf("---- Expected P_neut(i), 0 < i < n (neutral SFS per site) \n");
        for (i = 0; i < p->pm->n - 1; i++)
        {
            printf("E[P_neut(%d)] = %10.10f\n", i + 1, p->expec_neut[i]);
        }
        printf("\n---- Expected P_sel(i), 0 < i < n (selected SFS per site) \n");
        for (i = 0; i < p->pm->n - 1; i++)
        {
            printf("E[P_sel(%d)] = %10.10f\n", i + 1, p->expec_sel[i]);
        }

        // calculate observed divergence counts
        for (i = 0; i < p->no_neut; i++)
        {
            div_neut += p->counts_neut[i].sfs[p->pm->n - 1];
            len_neut += p->counts_neut[i].len_div;
        }
        for (i = 0; i < p->no_sel; i++)
        {
            div_sel += p->counts_sel[i].sfs[p->pm->n - 1];
            len_sel += p->counts_sel[i].len_div;
        }

        if (p->pm->div_flag == TRUE)
        {
            // print expected divergence counts
            printf("\n---- Expected D_neut and D_sel "
                   "(neutral and selected divergence per site) \n");
            printf("E[D_neut] = %10.10f\nE[D_sel] = %10.10f\n",
                   p->expec_neut[p->pm->n - 1], p->expec_sel[p->pm->n - 1]);
        }

        div_flag = p->pm->div_flag;

        // for the rest, I need p->pm->div_flag = TRUE
        p->pm->div_flag = TRUE;

        // in preparation for alpha, calculate the integral over the del DFE from eq. 8
        // and the integral over the full DFE from eq. 19
        // I can obtain this by getting the expected divergence with lambda = 1
        // and then later on subtract the misattributed polymorphism
        curr_lambda = p->pm->lambda;
        p->pm->lambda = 1;
        set_sel_expec(p->pm, &p->expec_sel, FALSE);
        expec_fixed_full = p->expec_sel[p->pm->n - 1];
        set_sel_expec(p->pm, &p->expec_sel, TRUE);
        expec_fixed_del = p->expec_sel[p->pm->n - 1];
        p->pm->lambda = curr_lambda;

        // calculate missatribuated polymorphism
        // I can do this easily by setting lambda to 0
        curr_lambda = p->pm->lambda;
        p->pm->lambda = 0;
        set_sel_expec(p->pm, &p->expec_sel, FALSE);
        mis_sel_full = p->expec_sel[p->pm->n - 1];
        set_sel_expec(p->pm, &p->expec_sel, TRUE);
        mis_sel_del = p->expec_sel[p->pm->n - 1];
        set_neut_expec(*p->pm, &p->expec_neut);
        mis_neut = p->expec_neut[p->pm->n - 1];

        // if I am not using divergence data, the below messages are pointless
        if (p->counts_neut[0].sfs[p->pm->n - 1] != -1)
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
        p->pm->lambda = curr_lambda;
        set_sel_lnL(p);
        set_neut_lnL(p);

        // remove misattributed polymorphism from expectations of fixed mutations
        expec_fixed_full -= mis_sel_full;
        expec_fixed_del -= mis_sel_del;

        // if I am not using divergence data, alpha_div should not be printed
        if (p->counts_neut[0].sfs[p->pm->n - 1] != -1)
        {
            // re-adjust divergence counts for the misattributed polymorhphism (eq. 9)
            div_neut -= len_neut * mis_neut;
            div_sel -= len_sel * mis_sel_full;
            // calculate alpha_div (eq. 8)
            printf("\n---- alpha_div = %.6f",
                   1 - len_sel * div_neut * expec_fixed_del / (len_neut * div_sel));
        }

        // calculate alpha_dfe (eq. 10)
        printf("\n---- alpha_dfe = %.6f\n\n",
                       1 - expec_fixed_del / expec_fixed_full);

        printf("\n\n");
    }

    // free memory
    gsl_vector_free(x);
    free(s);

    return (status);
}
