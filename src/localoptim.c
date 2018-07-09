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
#include <string.h>

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
    po->max_iter = 2000;
    po->minutes = 60;
    po->grad_descent = TRUE;
}

/****************************************************************************
 * Functions for printing
 ****************************************************************************/
void print_header_solution_neut(ParamsModel pm, int i, int with_fixed, FILE *f, int flag)
{
    char extra[10] = "";
    if (i >= 0)
    {
        sprintf(extra, "%d ", i);
    }
    char str[20];

    if (pm.eps_an_flag == flag
    		|| (with_fixed == TRUE && pm.eps_an_flag != SHARED))
    {
        strcpy(str, extra);
        strcat(str, "eps_an");
        fprintf(f, "%13.13s ", str);
    }
    if ((pm.lambda_flag == flag
    		|| (with_fixed == TRUE && pm.lambda_flag != SHARED)) && pm.div_flag == TRUE)
    {
        strcpy(str, extra);
        strcat(str, "lambda");
        fprintf(f, "%13.13s ", str);
    }
    if (pm.theta_bar_flag == flag
    		|| (with_fixed == TRUE && pm.theta_bar_flag != SHARED))
    {
        strcpy(str, extra);
        strcat(str, "theta_bar");
        fprintf(f, "%13.13s ", str);
    }
    if (pm.a_flag == flag
    		|| (with_fixed == TRUE && pm.a_flag != SHARED))
    {
        strcpy(str, extra);
        strcat(str, "a");
        fprintf(f, "%13.13s ", str);
    }
}

void print_header_solution_demo(ParamsModel pm, int i, int with_fixed, FILE *f, int flag)
{
    char extra[10] = "";
    if (i >= 0)
    {
        sprintf(extra, "%d ", i);
    }
    char str[20];

    unsigned j = 0;
    if (pm.r_flag == flag
    		|| (with_fixed == TRUE && pm.r_flag != SHARED))
    {
    	// print the group for each r parameter
    	int curr_r = 1;
    	for (j = 2; j <= pm.no_r; j++)
    	{
    		if (pm.inv_groups[curr_r] != pm.inv_groups[j])
    		{
    			strcpy(str, extra);
    			if (curr_r + 1 != j)
    			{
    				sprintf(str, "%s%s %d-%d", extra, "r", curr_r + 1, j);
    			}
    			else
    			{
    				sprintf(str, "%s%s %d", extra, "r", curr_r + 1);
    			}
    			fprintf(f, "%13.13s ", str);
    			curr_r = j;
    		}
    	}
//        for (j = 1; j < pm.no_groups; j++)
//        {
//            strcpy(str, extra);
//            sprintf(str, "%s%s%3d", extra, "r", j + 1);
//            fprintf(f, "%13.13s ", str);
//        }
    }
}

void print_header_solution_sel(ParamsModel pm, int i, int with_fixed, FILE *f, int flag)
{
    char extra[10] = "";
    if (i >= 0)
    {
        sprintf(extra, "%d ", i);
    }
    char str[20];
    unsigned j = 0;

    switch (pm.model)
    {
        case 1:
        {
            if (pm.sel_flag[0] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[0] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_bar");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[1] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[1] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "b");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[2] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[2] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_max");
                fprintf(f, "%13.13s ", str);
            }
            break;
        }
        case 2:
        {
            if (pm.sel_flag[0] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[0] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_d");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[1] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[1] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "b");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[2] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[2] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "p_b");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[3] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[3] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_b");
                fprintf(f, "%13.13s ", str);
            }
            break;
        }
        case 3:
        {
            if (pm.sel_flag[0] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[0] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_d");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[1] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[1] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "b");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[2] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[2] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "p_b");
                fprintf(f, "%13.13s ", str);
            }
            if (pm.sel_flag[3] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[3] != SHARED))
            {
                strcpy(str, extra);
                strcat(str, "S_b");
                fprintf(f, "%13.13s ", str);
            }
            break;
        }
        case 4:
        {
            for (j = 0; j < pm.no_sel/2; j++)
            {
                if (pm.sel_flag[2*j+1] == flag
                		|| (with_fixed == TRUE && pm.sel_flag[2*j+1] != SHARED))
                {
                    strcpy(str, extra);
                    sprintf(str, "%s%s", extra, "S_p");
                    if ((int) pm.sel_params[2*j] == pm.sel_params[2*j])
                    {
                        sprintf(str, "%s%4d", str, (int) pm.sel_params[2*j]);
                    }
                    else
                    {
                        sprintf(str, "%s%2.2g", str, pm.sel_params[2*j]);
                    }
                    fprintf(f, "%13.13s ", str);
                }
            }
            break;
        }
    }
}

void print_header_solution(ParamsModel pm, int i, int with_fixed, int flag,
                           int neut_ln, int sel_ln)
{
    if (neut_ln == TRUE)
    {
        print_header_solution_neut(pm, i, with_fixed, stdout, flag);
    }
    if (sel_ln == TRUE)
    {
            print_header_solution_sel(pm, i, with_fixed, stdout, flag);
    }
    if (neut_ln == TRUE)
    {
        print_header_solution_demo(pm, i, with_fixed, stdout, flag);
    }
}

void print_header_result(ParamsShare ps, int with_fixed, int grad_descent)
{
    unsigned i = 0;

    // I only print the iter for SHARED
    // and the lk, grad, status for TRUE
	printf("%4s ", "it");
	// if I am running on just one file, do not use i
	if (ps.no_data == 1)
	{
		print_header_solution(*ps.p[0].pm, -1, with_fixed, SHARED,
							  ps.use_neut_ln, ps.use_sel_ln);
	}
	else
	{
		print_header_solution(*ps.p[ps.which_r].pm, 0, with_fixed, SHARED,
							  ps.use_neut_ln, ps.use_sel_ln);
	}

	// if I am running on just one file, do not use i
	if (ps.no_data == 1)
	{
		print_header_solution(*ps.p[0].pm, -1, with_fixed, TRUE,
							  ps.use_neut_ln, ps.use_sel_ln);
	}
	else
	{
		if (ps.which == -1)
		{
			for (i = 0; i < ps.no_data; i++)
					{
						print_header_solution(*ps.p[i].pm, i + 1, with_fixed, TRUE,
											  ps.use_neut_ln, ps.use_sel_ln);
					}
		}
		else
		{
			print_header_solution(*ps.p[ps.which].pm, ps.which + 1, with_fixed, TRUE,
								  ps.use_neut_ln, ps.use_sel_ln);
		}
	}
	printf("%15.15s ",  "ln lk");
	if (grad_descent == TRUE)
	{
		printf("%13.13s ", "grad");
	}
	else
	{
		printf("%13.13s ", "size");
	}
	printf("%4s \n", "status");
}

char *removeTrailingZeros(char *s)
{
    int len = strlen(s);

    // remove trailling zeros
    // only if s contains '.'
    char *p = strchr (s,'.');
    if (p != NULL)
    {
        while (len > 0 && s[len - 1] == '0')
        {
            len--;
            s[len] = '\0';
        }
        // if all decimals were zeros, remove "."
        if (s[len - 1] == '.')
        {
            s[len - 1] = '\0';
        }
    }

    return (s);
}

void print_with_space(char **s, double x, FILE *f, int small)
{
    if ((int) x == x)
    {
        sprintf((*s), "%d", (int) x);
    }
    else
    {
        sprintf((*s), "%13.13f", x);
    }

    if (small == TRUE)
    {
        fprintf(f, "%4s ", *s);
    }
    else if (small == FALSE)
    {
        fprintf(f, "%13.13s ", removeTrailingZeros(*s));
    }
    else
    {
        // large space for likelihood
        fprintf(f, "%15.15s ", removeTrailingZeros(*s));
    }
}

void print_solution_neut(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag)
{
    if (pm.eps_an_flag == flag
    		|| (with_fixed == TRUE && pm.eps_an_flag != SHARED))
    {
        print_with_space(&s, pm.eps_an, f, FALSE);
    }
    if ((pm.lambda_flag == flag
    		|| (with_fixed == TRUE && pm.lambda_flag != SHARED)) && pm.div_flag == TRUE)
    {
        print_with_space(&s, pm.lambda, f, FALSE);
    }
    if (pm.theta_bar_flag == flag
    		|| (with_fixed == TRUE && pm.theta_bar_flag != SHARED))
    {
        print_with_space(&s, pm.theta_bar, f, FALSE);
    }
    if (pm.a_flag == flag
    		|| (with_fixed == TRUE && pm.a_flag != SHARED))
    {
        print_with_space(&s, pm.a, f, FALSE);
    }
}

void print_solution_demo(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag)
{
    unsigned i = 0;

    if (pm.r_flag == flag
    		|| (with_fixed == TRUE && pm.r_flag != SHARED))
    {
        // r[0] should always be 1
        for (i = 1; i < pm.no_groups; i++)
        {
            print_with_space(&s, pm.r[i], f, FALSE);
        }
    }
}

void print_solution_sel(ParamsModel pm, int with_fixed, char *s, FILE *f, int flag)
{
    unsigned i = 0;

    if (pm.model < 4)
    {
        for (i = 0; i < pm.no_sel; i++)
        {
            if (pm.sel_flag[i] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[i] != SHARED))
            {
                print_with_space(&s, pm.sel_params[i], f, FALSE);
            }
        }
    }
    else
    {
        for (i = 0; i < pm.no_sel/2; i++)
        {
            if (pm.sel_flag[2*i+1] == flag
            		|| (with_fixed == TRUE && pm.sel_flag[2*i+1] != SHARED))
            {
                print_with_space(&s, pm.sel_params[2*i+1], f, FALSE);
            }
        }
    }
}

void print_solution(ParamsModel pm, int with_fixed, char *s, int flag,
                    int neut_ln, int sel_ln)
{
    if (neut_ln == TRUE)
    {
            print_solution_neut(pm, with_fixed, s, stdout, flag);
    }
    if (sel_ln == TRUE)
    {
            print_solution_sel(pm, with_fixed, s, stdout, flag);
    }
    if (neut_ln == TRUE)
    {
        print_solution_demo(pm, with_fixed, s, stdout, flag);
    }
}

void print_result(ParamsShare ps, int with_fixed, int iter,
		          int status, char *s, double f, double grad)
{
    unsigned i = 0;

    // I only print the iter for SHARED
    // and the lk, grad, status for TRUE
    print_with_space(&s, iter, stdout, TRUE);
	print_solution(*ps.p[ps.which_r].pm, with_fixed, s, SHARED,
				   ps.use_neut_ln, ps.use_sel_ln);

	if (ps.which == -1)
	{
		for (i = 0; i < ps.no_data; i++)
		{
			print_solution(*ps.p[i].pm, with_fixed, s, TRUE,
					ps.use_neut_ln, ps.use_sel_ln);
		}
	}
	else
	{
		print_solution(*ps.p[ps.which].pm, with_fixed, s, TRUE,
					   ps.use_neut_ln, ps.use_sel_ln);
	}

	// print lk
	if (is_nan(f) == FALSE)
	{
		print_with_space(&s, -f, stdout, -200);
	}
	else
	{
		printf("%15.15s ", "NAN");
	}

	// print grad
	if (is_nan(grad) == FALSE)
	{
		print_with_space(&s, grad, stdout, FALSE);
	}
	else
	{
		printf("%13.13s ", "NAN");
	}

	print_with_space(&s, status, stdout, TRUE);
	printf("\n");
}

void print_result_gen(ParamsShare ps, int with_fixed, int iter, void *state,
		          int status, char *s, int grad_descent)
{
	// call print_result with the right arguments
	// according to the type of optimization
	double f = 0;
	double grad = -1;
	if (grad_descent == TRUE)
	{
		f = ((gsl_multimin_fdfminimizer *) state)->f;
		grad = gsl_blas_dnrm2(((gsl_multimin_fdfminimizer *) state)->gradient);
		// make sure ps contains the correct values
		undo_transform(&ps, ((gsl_multimin_fdfminimizer *) state)->x);
	}
	else
	{
		f = ((gsl_multimin_fminimizer *) state)->fval;
		// the gradient now contains the size
		grad = gsl_multimin_fminimizer_size(state);
		// make sure ps contains the correct values
		undo_transform(&ps, ((gsl_multimin_fminimizer *) state)->x);
	}
	print_result(ps, with_fixed, iter, status, s, f, grad);
}

/****************************************************************************
 * Functions for GSL optimization
 ****************************************************************************/
int set_lnL_f(const gsl_vector *x, void *pv, gsl_vector *f)
{
    gsl_vector_set(f, 0, get_lnL_share(x, pv));
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
    *f = get_lnL_share(x, pv);
    set_lnL_df(x, pv, df);
}

int optimize_partial_ln(const void *type,
                        ParamsOptim po, ParamsShare *ps, char *s)
{
	// I optimize either the neutral of the selected likelihood
    unsigned iter = 0;
    unsigned it = 0;

    // keep time of the running time
    // and stop when I exceed the allowed minutes
    clock_t begin = 0, end = 0;
    double time_spent = 0;
    int minutes = 0;
    begin = clock();

    int status = GSL_CONTINUE;
    int restart = FALSE;
    int same_state = 0;

    int size = 0;
    if (ps->use_neut_ln == TRUE)
    {
    	if (ps->which == -1)
    	{
    		size += ps->neut;
    	}
    	else
    	{
    		size += ps->p[ps->which].pm->neut;
    	}
    }
    if (ps->use_sel_ln == TRUE)
    {
    	if (ps->which == -1)
		{
			size += ps->sel;
		}
		else
		{
			size += ps->p[ps->which].pm->sel;
		}
    }

    // before initialization, make sure the counter is 0!
    ps->counter = 0;

    gsl_vector *x = gsl_vector_alloc(size);
    gsl_vector *prev_x = gsl_vector_alloc(size);
    // store best solution found before restarting - to avoid looping
    gsl_vector *prev_restart_x = gsl_vector_alloc(size);

    // initialize x
    transform(&x, *ps);
    gsl_vector_memcpy(prev_x, x);
    gsl_vector_memcpy(prev_restart_x, x);

    // depending on the type of optimization,
    // I have to use different functions from GSL
    void *state = NULL;
    void *func = NULL;
    gsl_vector *step_size = NULL;
    if (po.grad_descent == TRUE)
    {
    	func = malloc(sizeof(gsl_multimin_function_fdf));
    	((gsl_multimin_function_fdf *) func)->n = x->size;
    	((gsl_multimin_function_fdf *) func)->f = &get_lnL_share;
    	((gsl_multimin_function_fdf *) func)->df = &set_lnL_df;
    	((gsl_multimin_function_fdf *) func)->fdf = &set_lnL_fdf;
    	((gsl_multimin_function_fdf *) func)->params = (void *) ps;

    	// initialize state
		state = gsl_multimin_fdfminimizer_alloc(type, x->size);
		gsl_multimin_fdfminimizer_set(state, func, x, po.step_size, po.tol);
    }
    else
    {
    	func = malloc(sizeof(gsl_multimin_function));
    	((gsl_multimin_function *) func)->n = x->size;
    	((gsl_multimin_function *) func)->f = &get_lnL_share;
    	((gsl_multimin_function *) func)->params = (void *) ps;

    	// initialize state
		state = gsl_multimin_fminimizer_alloc(type, x->size);
		gsl_vector *step_size = gsl_vector_alloc(x->size);
		gsl_vector_set_all(step_size, po.step_size);
		gsl_multimin_fminimizer_set(state, func, x, step_size);
		// set the initial value of the function
		((gsl_multimin_fminimizer *) state)->fval = get_lnL_share(x, ps);
    }

    print_header_result(*ps, FALSE, po.grad_descent);
    print_result_gen(*ps, FALSE, 0, state, status, s, po.grad_descent);

    // do not start if lk is nan
    if (po.grad_descent == TRUE)
    {
    	if (is_nan(((gsl_multimin_fdfminimizer *) state)->f) == TRUE)
    	{
    		status = FOUND_NAN;
    	}
    }
    else
    {
    	if (is_nan(((gsl_multimin_fminimizer *) state)->fval) == TRUE)
    	{
    		status = FOUND_NAN;
    	}
    }

    // TODO
    // set max counter to 1000
    ps->max_counter = 1000;
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
        if (po.minutes > 0 && minutes > po.minutes)
        {
            status = OUT_OF_TIME;
            break;
        }

        // do optimization iteration
        // first, reset counter
        ps->counter = 0;

        if (po.grad_descent == TRUE)
        {
        	status = gsl_multimin_fdfminimizer_iterate(state);
        }
        else
        {
        	status = gsl_multimin_fminimizer_iterate(state);
        }

        if (status == GSL_EBADFUNC)
        {
            break;
        }

        // check the counter
        if (ps->counter >= ps->max_counter)
        {
        	status = MAX_LK;
            restart = TRUE;
        }
        else
        {
            // check if lk or gradient is nan
        	if (po.grad_descent == TRUE)
        	{
        		ps->criteria = gsl_blas_dnrm2(((gsl_multimin_fdfminimizer *) state)->gradient);
        		if (is_nan(((gsl_multimin_fdfminimizer *) state)->f) == TRUE
        				|| is_nan(ps->criteria) == TRUE)
        		{
        			status = FOUND_NAN;
        			break;
        		}
        	}
        	else
        	{
        		if (is_nan(((gsl_multimin_fminimizer *) state)->fval) == TRUE)
        		{
        			status = FOUND_NAN;
        			break;
        		}
        	}
        }

        // copy current solution to x
        if (po.grad_descent == TRUE)
		{
			gsl_vector_memcpy(x, ((gsl_multimin_fdfminimizer *) state)->x);
		}
		else
		{
			gsl_vector_memcpy(x, ((gsl_multimin_fminimizer *) state)->x);
		}

        // re-set max counter to 1000
        ps->max_counter = 1000;

        if (status)
        {
            // I am trying to restart
            // check if current solution
            // is different than previously found solution
            // but it took at least one optimizing cycle
            if (gsl_vector_equal(prev_restart_x, x) && it != 1 && status != MAX_LK)
            {
                // they're equal, stop optimization
                status = SAME_STATE_CIR;
                break;
            }
            else
            {
                gsl_vector_memcpy(prev_restart_x, x);
            }
            restart = TRUE;
        }

        // test for convergence
        if (status != MAX_LK)
        {
        	if (po.grad_descent == TRUE)
        	{
        		status = gsl_multimin_test_gradient(
							((gsl_multimin_fdfminimizer *) state)->gradient,
							po.eps_abs);
        	}
        	else
        	{
        		status = gsl_multimin_test_size(
        				gsl_multimin_fminimizer_size(state), po.eps_abs);
        	}
            if (status == GSL_SUCCESS)
            {
                break;
            }
        }

        // test if new solution is different than previous solution
        // for simplex, they can be the same for some reason,
        // but the algorithm still progresses
        if (po.grad_descent == TRUE && gsl_vector_equal(prev_x, x))
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
            gsl_vector_memcpy(prev_x, x);
            same_state = 0;
        }

        if (restart == TRUE)
        {
            if (status != MAX_LK)
            {
                status = RESTART;
            }
        }

        if (ps->verbose > 0 && iter % ps->verbose == 0)
        {
        	// if I had MAX_LK, I want to re-calculate the lk
        	if (status == MAX_LK)
        	{
            	// re-evaluate the likelihood
                ps->counter = 0;
                if (po.grad_descent == TRUE)
                {
                	((gsl_multimin_fdfminimizer *) state)->f = get_lnL_share(x, ps);
                }
                else
                {
                	((gsl_multimin_fminimizer *) state)->fval = get_lnL_share(x, ps);
                }

        	}
            print_result_gen(*ps, FALSE, iter, state, status, s, po.grad_descent);
        }

        if (restart == TRUE)
        {
            ps->counter = 0;
            // restart the optimization
            // the gsl restart function does not seem to work properly
            // free state and reallocate it
            if (po.grad_descent)
            {
            	gsl_multimin_fdfminimizer_free(state);
				state = gsl_multimin_fdfminimizer_alloc(type, x->size);
				gsl_multimin_fdfminimizer_set(state, func, x, po.step_size, po.tol);
            }
            else
            {
            	gsl_multimin_fminimizer_free(state);
            	state = gsl_multimin_fminimizer_alloc(type, x->size);
				step_size = gsl_vector_alloc(x->size);
				gsl_vector_set_all(step_size, po.step_size);
				gsl_multimin_fminimizer_set(state, func, x, step_size);
				// set the initial value of the function
				((gsl_multimin_fminimizer *) state)->fval = get_lnL_share(x, ps);
            }

            // if I am restarting because I reached MAX_LK
            // then I increase the max counter
            if (status == MAX_LK)
            {
                ps->max_counter = 200 * size;
                // though never set it to less than 1000
                ps->max_counter = ps->max_counter < 1000 ? 1000 : ps->max_counter;
                ps->max_counter = 1000;
            }
            restart = FALSE;
            it = 0;
        }
    }

    // if I had MAX_LK, I want to re-calculate the lk and grad
	if (status == MAX_LK)
	{
		ps->counter = 0;
		ps->max_counter = 1000;
		// free state and reallocate it
		if (po.grad_descent)
		{
			gsl_multimin_fdfminimizer_free(state);
			state = gsl_multimin_fdfminimizer_alloc(type, x->size);
			gsl_multimin_fdfminimizer_set(state, func, x, po.step_size, po.tol);
		}
		else
		{
			gsl_multimin_fminimizer_free(state);
			state = gsl_multimin_fminimizer_alloc(type, x->size);
			step_size = gsl_vector_alloc(x->size);
			gsl_vector_set_all(step_size, po.step_size);
			gsl_multimin_fminimizer_set(state, func, x, step_size);
		}
	}

	// before printing final result
    // make sure I store best solution found
    if (po.grad_descent == TRUE)
	{
		gsl_vector_memcpy(x, gsl_multimin_fdfminimizer_x(state));
		gsl_vector_memcpy(((gsl_multimin_fdfminimizer *) state)->x, x);
	    ((gsl_multimin_fdfminimizer *) state)->f = gsl_multimin_fdfminimizer_minimum(state);
	    gsl_vector_memcpy(((gsl_multimin_fdfminimizer *) state)->gradient,
	    		gsl_multimin_fdfminimizer_gradient(state));
	    ps->criteria = gsl_blas_dnrm2(((gsl_multimin_fdfminimizer *) state)->gradient);
	}
	else
	{
		gsl_vector_memcpy(x, gsl_multimin_fminimizer_x(state));
		gsl_vector_memcpy(((gsl_multimin_fminimizer *) state)->x, x);
		((gsl_multimin_fminimizer *) state)->fval = gsl_multimin_fminimizer_minimum(state);
		ps->criteria = gsl_multimin_fminimizer_size(state);
	}
    get_lnL_share(x, ps);
    print_result_gen(*ps, FALSE, iter, state, status, s, po.grad_descent);

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
        printf("-- Likelihood or gradient is NAN\n");
    } else if (status == OUT_OF_TIME)
    {
        printf("-- Local optimization: reached maximum running time allowed\n");
    }
    else if (status == MAX_LK)
    {
        printf("-- Local optimization: reached maximum number of likelihood "
        		"evaluations allowed.\n");
    }
    else if (status)
    {
        printf("-- Local optimization: %s\n", gsl_strerror(status));
    }

    // free memory
    free(func);
    gsl_vector_free(x);
    gsl_vector_free(prev_x);
    gsl_vector_free(prev_restart_x);
    if (po.grad_descent == TRUE)
    {
    	gsl_multimin_fdfminimizer_free(state);
    }
    else
    {
    	gsl_multimin_fminimizer_free(state);
    }

    // make sure the counter is reset to 0
    ps->counter = 0;
    ps->max_counter = 1000;

    return (EXIT_SUCCESS);
}

int optimize_per_file(const void *type,
		              ParamsOptim po, ParamsShare *ps, char *s)
{
	// if no parameters are shared between the files
	// I need to run the optimization for each file separately

	// the easy one, I have shared parameters or just one file
	ps->which = -1;
	if (ps->no_data == 1 || ps->idx[0] > 0)
	{
		optimize_partial_ln(type, po, ps, s);
		return (EXIT_SUCCESS);
	}

	// run each file separately
	size_t i = 0;
	for (i = 0; i < ps->no_data; i++)
	{
		ps->which = i;
		optimize_partial_ln(type, po, ps, s);

	}

	// I need to calculate the joint likelihood and gradient
	ps->which = -1;
	int size = 0;
	if (ps->use_neut_ln == TRUE)
	{
		size += ps->neut;
	}
	if (ps->use_sel_ln == TRUE)
	{
		size += ps->sel;
	}

	gsl_vector *x = gsl_vector_alloc(size);
	// initialize x
	transform(&x, *ps);

	// before initialization, make sure the counter is 0!
	ps->counter = 0;

	// calculate joint lk and gradient / size
	if (po.grad_descent == TRUE)
	{
		gsl_multimin_function_fdf func;
		func.n = x->size;
		func.f = &get_lnL_share;
		func.df = &set_lnL_df;
		func.fdf = &set_lnL_fdf;
		func.params = (void *) ps;

		// initialize state
		gsl_multimin_fdfminimizer *state = gsl_multimin_fdfminimizer_alloc(type,
				                                                       x->size);
		gsl_multimin_fdfminimizer_set(state, &func, x, po.step_size, po.tol);
		ps->criteria = gsl_blas_dnrm2(state->gradient);

		// print the info on the joint likelihood and gradient
		if (ps->use_neut_ln == TRUE && ps->use_sel_ln == TRUE)
		{
			printf("-- Joint likelihood over all files %.15f and gradient %.5f\n",
					-state->f, ps->criteria);
		}
		else
		{
			if (ps->use_neut_ln == TRUE)
			{
				printf("-- Joint neutral likelihood over all files %.15f and gradient %.5f\n",
						-state->f, ps->criteria);
			}
			if (ps->use_sel_ln == TRUE)
			{
				printf("-- Joint selected likelihood over all files %.15f and gradient %.5f\n",
						-state->f, ps->criteria);
			}
		}

		// free memory
		gsl_multimin_fdfminimizer_free(state);
	}
	else
	{
		gsl_multimin_function func;
		func.n = x->size;
		func.f = &get_lnL_share;
		func.params = (void *) ps;

		// initialize state
		gsl_multimin_fminimizer *state = gsl_multimin_fminimizer_alloc(type,
				                                                       x->size);
		gsl_vector *step_size = gsl_vector_alloc(x->size);
		gsl_vector_set_all(step_size, po.step_size);
		gsl_multimin_fminimizer_set(state, &func, x, step_size);
		ps->criteria = gsl_multimin_fminimizer_size(state);
		// set the initial value of the function
		state->fval = get_lnL_share(x, ps);

		// print the info on the joint likelihood and gradient
		if (ps->use_neut_ln == TRUE && ps->use_sel_ln == TRUE)
		{
			printf("-- Joint likelihood over all files %.15f and size %.5f\n",
					-state->fval, ps->criteria);
		}
		else
		{
			if (ps->use_neut_ln == TRUE)
			{
				printf("-- Joint neutral likelihood over all files %.15f and size %.5f\n",
						-state->fval, ps->criteria);
			}
			if (ps->use_sel_ln == TRUE)
			{
				printf("-- Joint selected likelihood over all files %.15f and size %.5f\n",
						-state->fval, ps->criteria);
			}
		}

		// free memory
		gsl_multimin_fminimizer_free(state);
	}

	// free memory
	gsl_vector_free(x);

	// store best optimum found
	ps->lnL = - ps->lnL_neut - ps->lnL_sel;

	return (EXIT_SUCCESS);
}

int optimize(const void *type, ParamsOptim po, ParamsShare *ps, char *s)
{
    // ps->kind: kind of likelihood to optimize
    // 0: neutral, selected
    // 1: neutral, selected, joint
    // 2: joint

    if (ps->kind < 2)
    {
        if (ps->neut > 0)
        {
            // first, optimize the neutral likelihood
            printf("-- Optimizing neutral parameters\n");
            ps->use_neut_ln = TRUE;
            ps->use_sel_ln = FALSE;
            optimize_per_file(type, po, ps, s);
        }
        if (ps->sel > 0)
        {
            // then optimize the selected likelihood
            printf("-- Optimizing selected parameters\n");
            ps->use_neut_ln = FALSE;
            ps->use_sel_ln = TRUE;
            optimize_per_file(type, po, ps, s);
        }
    }

    if (ps->kind == 0)
    {
    	// I need to calculate the joint gradient / size
    	ps->use_neut_ln = TRUE;
    	ps->use_sel_ln = TRUE;
    	int size = ps->neut + ps->sel;

    	gsl_vector *x = gsl_vector_alloc(size);
    	// initialize x
    	transform(&x, *ps);

    	// calculate joint lk and gradient / size
    	if (po.grad_descent == TRUE)
    	{
    		gsl_multimin_function_fdf func;
    		func.n = x->size;
    		func.f = &get_lnL_share;
    		func.df = &set_lnL_df;
    		func.fdf = &set_lnL_fdf;
    		func.params = (void *) ps;

    		// initialize state
    		gsl_multimin_fdfminimizer *state = gsl_multimin_fdfminimizer_alloc(type,
    				x->size);
    		gsl_multimin_fdfminimizer_set(state, &func, x, po.step_size, po.tol);
    		ps->criteria = gsl_blas_dnrm2(state->gradient);

    		// print the info on the joint likelihood and gradient
    		printf("-- Joint neutral and selected likelihood %.15f and gradient %.5f\n",
    				-state->f, ps->criteria);

    		// free memory
    		gsl_multimin_fdfminimizer_free(state);
    	}
    	else
    	{
    		gsl_multimin_function func;
    		func.n = x->size;
    		func.f = &get_lnL_share;
    		func.params = (void *) ps;

    		// initialize state
    		gsl_multimin_fminimizer *state = gsl_multimin_fminimizer_alloc(type,
    				x->size);
    		gsl_vector *step_size = gsl_vector_alloc(x->size);
    		gsl_vector_set_all(step_size, po.step_size);
    		gsl_multimin_fminimizer_set(state, &func, x, step_size);
    		ps->criteria = gsl_multimin_fminimizer_size(state);
    		// set the initial value of the function
    		state->fval = get_lnL_share(x, ps);

    		// print the info on the joint likelihood and gradient
    		printf("-- Joint neutral and selected likelihood %.15f and size %.5f\n",
    				-state->fval, ps->criteria);

    		// free memory
    		gsl_multimin_fminimizer_free(state);
    	}

    	// free memory
    	gsl_vector_free(x);
    }

    if (ps->kind > 0)
    {
        // now optimize jointly
        printf("-- Optimizing all parameters\n");
        ps->use_neut_ln = TRUE;
        ps->use_sel_ln = TRUE;
        optimize_per_file(type, po, ps, s);
    }

    // store best optimum found
    ps->lnL = - ps->lnL_neut - ps->lnL_sel;

    return (EXIT_SUCCESS);
}
