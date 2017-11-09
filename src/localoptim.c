/*
 * polyDFE v1.0: predicting DFE and alpha from polymorphism data
 * Copyright (c) 2016  Paula Tataru and Marco A.P. Franco
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

#include "localoptim.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_nan.h>
#include <gsl/gsl_sys.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "transform.h"

int is_nan(double x)
{
    if (gsl_isnan(x) == 1 || x == DBL_MAX)
    {
        return (TRUE);
    }
    return (FALSE);
}

/****************************************************************************
 * Functions on structures: initialization
 ****************************************************************************/
void initialize_params_optim(ParamsOptim *po)
{
    po->eps_abs = 0.0001;
    po->step_size = 2;
    po->tol = 0.1;
    po->max_iter = 1500;
    po->minutes = 60;
}

/****************************************************************************
 * Functions for printing
 ****************************************************************************/
void print_hearder_solution_neut(ParamsModel pm, int with_fixed, FILE *f)
{
    if (pm.eps_an_flag == TRUE || with_fixed == TRUE)
    {
        fprintf(f, "     eps_an   ");
    }
    if (pm.eps_cont_flag == TRUE || with_fixed == TRUE)
    {
        fprintf(f, "     eps_cont ");
    }
    if ((pm.lambda_flag == TRUE || with_fixed == TRUE) && pm.div_flag == TRUE)
    {
        fprintf(f, "     lambda   ");
    }
    if (pm.theta_bar_flag == TRUE || with_fixed == TRUE)
    {
        fprintf(f, "     theta_bar");
    }
    if (pm.a_flag == TRUE || with_fixed == TRUE)
    {
        fprintf(f, "     a        ");
    }
}

void print_hearder_solution_demo(ParamsModel pm, int with_fixed, FILE *f)
{
    unsigned i = 0;
    if (pm.r_flag == TRUE || with_fixed == TRUE)
    {
        for (i = 1; i < pm.no_groups; i++)
        {
            fprintf(f, "     r%3d     ", i + 1);
        }
    }
}

void print_hearder_solution_sel(ParamsModel pm, int with_fixed, FILE *f)
{
    unsigned i = 0;

    switch (pm.model)
    {
        case 1:
        {
            if (pm.sel_flag[0] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "      S_bar     ");
            }
            if (pm.sel_flag[1] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     b        ");
            }
            if (pm.sel_flag[2] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     S_max    ");
            }
            break;
        }
        case 2:
        {
            if (pm.sel_flag[0] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "      S_d       ");
            }
            if (pm.sel_flag[1] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     b        ");
            }
            if (pm.sel_flag[2] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     p_b      ");
            }
            if (pm.sel_flag[3] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     S_b      ");
            }
            break;
        }
        case 3:
        {
            if (pm.sel_flag[0] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "      S_d       ");
            }
            if (pm.sel_flag[1] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     b        ");
            }
            if (pm.sel_flag[2] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     p_b      ");
            }
            if (pm.sel_flag[3] == TRUE || with_fixed == TRUE)
            {
                fprintf(f, "     S_b      ");
            }
            break;
        }
        case 4:
        {
            for (i = 0; i < pm.no_sel/2; i++)
            {
                if (pm.sel_flag[2*i+1] != FALSE || with_fixed == TRUE)
                {
                    if ((int) pm.sel_params[2*i] == pm.sel_params[2*i])
                    {
                        fprintf(f, "     S_p%4d  ", (int) pm.sel_params[2*i]);
                    }
                    else
                    {
                        fprintf(f, "     S_p%2.2g", pm.sel_params[2*i]);
                    }
                }
            }
            break;
        }
    }
}

void print_hearder_solution(ParamsModel pm, int with_fixed)
{
    print_hearder_solution_neut(pm, with_fixed, stdout);
    print_hearder_solution_sel(pm, with_fixed, stdout);
    print_hearder_solution_demo(pm, with_fixed, stdout);
}

void print_header_result(ParamsModel pm, int with_fixed)
{
    printf(" it  ");
    print_hearder_solution(pm, with_fixed);
    printf("    ln lk             grad     status\n");
}

void print_with_space(char **s, double x, FILE *f)
{
    char *aux;

    if (x < 10 && x > -10)
    {
        aux = "     ";
    }
    else if (x < 100 && x > -100)
    {
        aux = "    ";
    }
    else if (x < 1000 && x > -1000)
    {
        aux = "   ";
    }
    else if (x < 10000 && x > -10000)
    {
        aux = "  ";
    }
    else
    {
        aux = " ";
    }

    if (x >= 0.00001 || -x >= 0.00001)
    {
        sprintf((*s), "%s%.5f", aux, x);
    } else
    {
        aux = "    ";
        sprintf((*s), "%s%.2e", aux, x);
    }

    fprintf(f, "%s  ", *s);
}

void print_solution_neut(ParamsModel pm, int with_fixed, gsl_vector *x, char *s,
                         FILE *f)
{
    if (x)
    {
        undo_transform(&pm, x);
    }

    if (pm.eps_an_flag == TRUE || with_fixed == TRUE)
    {
        print_with_space(&s, pm.eps_an, f);
    }
    if (pm.eps_cont_flag == TRUE || with_fixed == TRUE)
    {
        print_with_space(&s, pm.eps_cont, f);
    }
    if ((pm.lambda_flag == TRUE || with_fixed == TRUE) && pm.div_flag == TRUE)
    {
        print_with_space(&s, pm.lambda, f);
    }
    if (pm.theta_bar_flag == TRUE || with_fixed == TRUE)
    {
        print_with_space(&s, pm.theta_bar, f);
    }
    if (pm.a_flag == TRUE || with_fixed == TRUE)
    {
        print_with_space(&s, pm.a, f);
    }
}

void print_solution_demo(ParamsModel pm, int with_fixed, gsl_vector *x, char *s,
                         FILE *f)
{
    unsigned i = 0;
    if (x)
    {
        undo_transform(&pm, x);
    }

    if (pm.r_flag == TRUE || with_fixed == TRUE)
    {
        // r[0] should always be 1
        for (i = 1; i < pm.no_groups; i++)
        {
            print_with_space(&s, pm.r[i], f);
        }
    }
}

void print_solution_sel(ParamsModel pm, int with_fixed, gsl_vector *x, char *s,
                        FILE *f)
{
    unsigned i = 0;
    if (x)
    {
        undo_transform(&pm, x);
    }

    if (pm.model < 4)
    {
        for (i = 0; i < pm.no_sel; i++)
        {
            if (pm.sel_flag[i] == TRUE || with_fixed == TRUE)
            {
                print_with_space(&s, pm.sel_params[i], f);
            }
        }
    }
    else
    {
        for (i = 0; i < pm.no_sel/2; i++)
        {
            if (pm.sel_flag[2*i+1] != FALSE || with_fixed == TRUE)
            {
                print_with_space(&s, pm.sel_params[2*i+1], f);
            }
        }
    }
}

void print_solution(ParamsModel pm, int with_fixed, gsl_vector *x, char *s)
{
    print_solution_neut(pm, with_fixed, x, s, stdout);
    print_solution_sel(pm, with_fixed, x, s, stdout);
    print_solution_demo(pm, with_fixed, x, s, stdout);
}

void print_result(ParamsModel pm, int with_fixed, int iter,
                  gsl_multimin_fdfminimizer *state, int status, char *s)
{
    double grad = gsl_blas_dnrm2(state->gradient);
    printf("%3d  ", iter);
    print_solution(pm, with_fixed, state->x, s);
    // print nan if necessary
    if (is_nan(state->f) == FALSE)
    {
        printf("%.10f ", -state->f);
        // print nan if necessary
		if (is_nan(grad) == FALSE)
		{
			print_with_space(&s, grad, stdout);
		}
		else
		{
			printf("   NAN ");
		}
    }
    else
    {
        printf("     NAN               NAN      ");
    }
	printf("%3d\n", status);
}

/****************************************************************************
 * Functions for GSL optimization
 ****************************************************************************/
int set_lnL_f(const gsl_vector *x, void *pv, gsl_vector *f)
{
    gsl_vector_set(f, 0, get_lnL(x, pv));
    return (GSL_SUCCESS);
}

void set_lnL_df(const gsl_vector *x, void *pv, gsl_vector *df)
{
    // need to calculate partial derivatives of the log likelihood
    // old code used the Jacobian matrix

    gsl_multifit_function_fdf fdf;
    fdf.f = &set_lnL_f;
    fdf.df = NULL;
    fdf.fdf = NULL;
    fdf.n = 1;
    fdf.p = x->size;
    fdf.params = pv;
    // GSL 2.1 version
    // fdf.nevalf = 1;
    // fdf.nevaldf = 1;

    gsl_matrix *J = gsl_matrix_alloc(1, fdf.p);
    J->size1 = fdf.n;
    J->size2 = fdf.p;

    gsl_vector *f = gsl_vector_alloc(1);
    set_lnL_f(x, pv, f);

    // GSL 2.1 version - but it behaves weirdly
    // int status = gsl_multifit_fdfsolver_dif_df(x, NULL, &fdf, f, J);
    // GSL 1.16 version
    int status = gsl_multifit_fdfsolver_dif_df(x, &fdf, f, J);

    if (status)
    {
        // gsl_vector_set_all(df, GSL_NAN);
        gsl_vector_set_all(df, DBL_MAX);
    }
    else
    {
        gsl_matrix_get_row(df, J, 0);
    }

    gsl_vector_free(f);
    gsl_matrix_free(J);
}

void set_lnL_fdf(const gsl_vector *x, void *pv, double *f, gsl_vector *df)
{
    *f = get_lnL(x, pv);
    set_lnL_df(x, pv, df);
}

int optimize_using_derivatives(const gsl_multimin_fdfminimizer_type *type,
                               ParamsOptim po, Params *p, char *s)
{
    unsigned iter = 0;
    unsigned it = 0;
    double grad = 0;

    // keep time of the running time
    // and stop when I exceed the allowed minutes
    clock_t begin = 0, end = 0;
    double time_spent = 0;
    int minutes = 0;
    begin = clock();

    int status = GSL_CONTINUE;
    int restart = FALSE;
    int same_state = 0;

    int size = p->pm->neut + p->pm->sel;
    gsl_vector *x = gsl_vector_alloc(size);
    gsl_vector *prev_x = gsl_vector_alloc(size);
    // store best solution found before restarting - to avoid looping
    gsl_vector *prev_restart_x = gsl_vector_alloc(size);

    // initialize x
    transform(&x, *p->pm);
    gsl_vector_memcpy(prev_x, x);
    gsl_vector_memcpy(prev_restart_x, x);

    gsl_multimin_function_fdf func;
    func.n = x->size;
    func.f = &get_lnL;
    func.df = &set_lnL_df;
    func.fdf = &set_lnL_fdf;
    func.params = (void *) p;

    // initialize state
    gsl_multimin_fdfminimizer *state = gsl_multimin_fdfminimizer_alloc(type,
                                                                       x->size);
    gsl_multimin_fdfminimizer_set(state, &func, x, po.step_size, po.tol);

    print_header_result(*p->pm, FALSE);
    print_result(*p->pm, FALSE, 0, state, status, s);

    // do not start if lk is nan
    if (is_nan(state->f) == TRUE)
    {
        status = FOUND_NAN;
    }
    for (iter = 1;
                    iter <= po.max_iter
                        && (status == GSL_CONTINUE
                                        || status == RESTART
                                        || status == MAX_LK);
                    iter++, it++)
    {
        // check the running time
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        minutes = floor(time_spent / 60);
        if (minutes > po.minutes)
        {
            status = OUT_OF_TIME;
            break;
        }

        // do optimization iteration
        // first, reset counter
        p->counter = 0;
        status = gsl_multimin_fdfminimizer_iterate(state);
        if (status == GSL_EBADFUNC)
        {
            break;
        }

        // check the counter
        if (p->counter >= p->max_counter)
        {
            status = MAX_LK;
            // restart the optimization
            restart = TRUE;
        }
        else
        {

            // check if lk or gradient is nan
            grad = gsl_blas_dnrm2(state->gradient);
            if (is_nan(state->f) == TRUE || is_nan(grad) == TRUE)
            {
                status = FOUND_NAN;
                break;
            }
        }

        if (status)
        {
            // I am trying to restart
            // check if current solution
            // is different than previously found solution
            // but it took at least one optimizing cycle
            if (gsl_vector_equal(prev_restart_x, state->x) && it != 1
                            && status != MAX_LK)
            {
                // they're equal, stop optimization
                status = SAME_STATE_CIR;
                break;
            }
            else
            {
                gsl_vector_memcpy(prev_restart_x, state->x);
            }
            restart = TRUE;
        }

        // test for convergence
        if (status != MAX_LK)
        {
            status = gsl_multimin_test_gradient(state->gradient, po.eps_abs);
            if (status == GSL_SUCCESS)
            {
                break;
            }
        }

        // test if new solution is different than previous solution
        if (gsl_vector_equal(prev_x, state->x))
        {
            restart = TRUE;
            same_state++;
            if (same_state == 3)
            {
                if (status != MAX_LK)
                {
                    status = SAME_STATE;
                }
                break;
            }
        }
        else
        {
            gsl_vector_memcpy(prev_x, state->x);
            same_state = 0;
        }

        if (restart == TRUE)
        {
            if (status != MAX_LK)
            {
                status = RESTART;
            }
        }

        if (p->verbose > 0 && iter % p->verbose == 0)
        {
            print_result(*p->pm, FALSE, iter, state, status, s);
        }

        if (restart == TRUE)
        {
            restart = FALSE;
            it = 0;
            // restart the optimization
            // the gsl restart function does not work on its own
            gsl_vector_memcpy(x, state->x);
            // recalculate likelihood
            state->f = get_lnL(x, p);
            gsl_multimin_fdfminimizer_set(state, &func, x,
                                          po.step_size,
                                          po.tol);
            gsl_multimin_fdfminimizer_restart(state);
        }
    }

    print_result(*p->pm, FALSE, iter, state, status, s);
    if (iter > po.max_iter)
    {
        printf("-- Local optimization: reached maximum number of iterations "
               "allowed\n");
    }
    else if (status == SAME_STATE)
    {
        printf("-- Local optimization: stuck in the same solution\n");
    }
    else if (status == SAME_STATE_CIR)
    {
        printf("-- Local optimization: restarting lead to the same solution\n");
    }
    else if (status == FOUND_NAN)
    {
        printf("-- Likelihood or gradient is nan\n");
    } else if (status == OUT_OF_TIME)
    {
        printf("-- Local optimization: reached maximum running time allowed\n");
    }
    else if (status == MAX_LK)
    {
        printf("-- Reached maximum number of likelihood evaluations allowed. "
               "Results are not reliable\n");
    }
    else if (status)
    {
        printf("-- Local optimization: %s\n", gsl_strerror(status));
    }

    // store best optimum found
    p->lnL = state->f;
    undo_transform(p->pm, state->x);
    p->grad = gsl_blas_dnrm2(state->gradient);

    // free memory
    gsl_vector_free(x);
    gsl_vector_free(prev_x);
    gsl_vector_free(prev_restart_x);
    gsl_multimin_fdfminimizer_free(state);

    return (EXIT_SUCCESS);
}
