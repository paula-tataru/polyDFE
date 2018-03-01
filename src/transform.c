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

#include "transform.h"

#include <gsl/gsl_vector_double.h>
#include <math.h>
#include <stdlib.h>

#include "likelihood.h"

/****************************************************************************
 * Generalized logistic functions
 * https://en.wikipedia.org/wiki/Generalised_logistic_function
 * to transform between y in (a, b) and x in (-inf, +inf)
 * when x = 0, y = y0
 * k determines how concentrated the values are around x0
 * smaller k means more spread
 * when y0 = (a+b)/2, spread is typically in interval [x0-5/k, x0+5/k]
 ****************************************************************************/
double get_logistic(double x, double a, double b, double y0, double k)
{
    double q = (b - a) / (y0 - a) - 1;
    return (a + (b - a) / (1 + q * exp(-k * x)));
}

double get_inverse_logistic(double y, double a, double b, double y0, double k)
{
    double q = (b - a) / (y0 - a) - 1;
    return (-log(((b - a) / (y - a) - 1) / q) / k);
}

/****************************************************************************
 * Transformation functions good for various ranges
 ****************************************************************************/
double tran(double x, double a, double b, double k)
{
    return (get_inverse_logistic(x, a, b, (a + b) / 2, k));
}

double untran(double y, double a, double b, double k)
{
    return (get_logistic(y, a, b, (a + b) / 2, k));
}

// for the step, I want to transform from (0, +inf) to (0, b) instead
// for that, I change the original transformation from (-inf, +inf)
double untran_step(double y, double b, double k)
{
    return (2 * (untran(y, 0, b, k) - b / 2));
}

// for the _one versions, x is positive, with a tipping point at 1, y0 = 1
double tran_one(double x, double a, double b, double k)
{
    // only transform around 1 if b > 1
    // otherwise use the standard transformation
    if (b > 1)
    {
        return (get_inverse_logistic(x, a, b, 1, k));
    }
    return (get_inverse_logistic(x, a, b, (a + b) / 2, k));
}

double untran_one(double y, double a, double b, double k)
{
    // only transform around 1 if b > 1
    // otherwise use the standard transformation
    if (b > 1)
    {
        return (get_logistic(y, a, b, 1, k));
    }
    return (get_logistic(y, a, b, (a + b) / 2, k));
}

/****************************************************************************
 * Transformation functions for model parameters
 ****************************************************************************/
void transform(gsl_vector **x, ParamsModel pm)
{
    unsigned i = 0;
    unsigned j = 0;

    if (pm.r_flag == TRUE)
    {
        // r[0] should always be 1
        for (i = 1; i < pm.no_groups; i++)
        {
            gsl_vector_set((*x), i - 1,
                           tran_one(pm.r[i], pm.r_min, pm.r_max, pm.k));
        }
        i--;
    }
    if (pm.lambda_flag == TRUE)
    {
        gsl_vector_set((*x), i++,
                       tran(pm.lambda, pm.lambda_min, pm.lambda_max, pm.k));
    }
    if (pm.theta_bar_flag == TRUE)
    {
        gsl_vector_set((*x),
                       i++,
                       tran(pm.theta_bar, pm.theta_bar_min, pm.theta_bar_max,
                            pm.k));
    }
    if (pm.a_flag == TRUE)
    {
        gsl_vector_set((*x),
                       i++,
                       tran_one(pm.a, pm.a_min, pm.a_max, pm.k));
    }
    if (pm.eps_an_flag == TRUE)
    {
        gsl_vector_set((*x), i++,
                       tran(pm.eps_an, pm.eps_an_min, pm.eps_an_max, pm.k));
    }
    if (pm.eps_cont_flag == TRUE)
    {
        gsl_vector_set((*x),
                       i++,
                       tran(pm.eps_cont, pm.eps_cont_min, pm.eps_cont_max,
                            pm.k));
    }


    // b - position 1 - should be transformed around one
    // everything else has the normal transformation
    // I don't need to treat model D separately because I have the flags set
    // to FALSE for the selection coefficients
    for (j = 0; j < pm.no_sel; j++)
    {
        if (pm.sel_flag[j] != TRUE)
        {
            // skip this
            continue;
        }
        if (pm.model < 4 && j == 1)
        {
            gsl_vector_set((*x), i++,
                           tran_one(pm.sel_params[j], pm.sel_min[j], pm.sel_max[j], pm.k));
        }
        else
        {
            gsl_vector_set((*x), i++,
                           tran(pm.sel_params[j], pm.sel_min[j], pm.sel_max[j], pm.k));
        }
    }

}

void undo_transform(ParamsModel *pm, const gsl_vector *x)
{
    unsigned i = 0;
    unsigned j = 0;

    if (pm->r_flag == TRUE)
    {
        // r[0] should always be 1
        for (i = 1; i < pm->no_groups; i++)
        {
            pm->r[i] = untran_one(gsl_vector_get(x, i - 1), pm->r_min,
                                  pm->r_max, pm->k);
        }
        i--;
    }
    if (pm->lambda_flag == TRUE)
    {
        pm->lambda = untran(gsl_vector_get(x, i++), pm->lambda_min,
                            pm->lambda_max, pm->k);
    }
    if (pm->theta_bar_flag == TRUE)
    {
        pm->theta_bar = untran(gsl_vector_get(x, i++), pm->theta_bar_min,
                               pm->theta_bar_max, pm->k);
    }
    if (pm->a_flag == TRUE)
    {
        pm->a = untran_one(gsl_vector_get(x, i++), pm->a_min,
                               pm->a_max, pm->k);
    }
    if (pm->eps_an_flag == TRUE)
    {
        pm->eps_an = untran(gsl_vector_get(x, i++), pm->eps_an_min,
                            pm->eps_an_max, pm->k);
    }
    if (pm->eps_cont_flag == TRUE)
    {
        pm->eps_cont = untran(gsl_vector_get(x, i++), pm->eps_cont_min,
                              pm->eps_cont_max, pm->k);
    }

    // b - position 1 - should be transformed around one
    // everything else has the normal transformation
    for (j = 0; j < pm->no_sel; j++)
    {
        if (pm->sel_flag[j] != TRUE)
        {
            // skip this
            continue;
        }
        if (pm->model < 4 && j == 1)
        {
            pm->sel_params[j] = untran_one(gsl_vector_get(x, i++), pm->sel_min[j], pm->sel_max[j],
                                           pm->k);
        }
        else
        {
            pm->sel_params[j] = untran(gsl_vector_get(x, i++), pm->sel_min[j], pm->sel_max[j],
                                       pm->k);
        }
    }
    // if model 4, set the first estimated parameter to 1 - sum(rest)
    if (pm->model == 4)
    {
        pm->sel_params[pm->sel_fixed] = 1;
        for (j = 0; j < pm->no_sel/2; j++)
        {
            if (2*j+1 != pm->sel_fixed)
            {
                pm->sel_params[pm->sel_fixed] -= pm->sel_params[2*j+1];
            }
        }
    }
}

double one_step(double step, double value, double min_value, double max_value)
{
    // randomly change value, in the range min_value and max_value,
    // using step
    // make sure step is positive - it can be negative in the case of s_bar
    step = step > 0 ? step : -step;
    // find a new value that is approximately +- step away from the current value
    double low = min_value >= value - step ? min_value : value - step;
    double up = max_value <= value + step ? max_value : value + step;
    double u = drand48();
    u = low + (up - low) * u;
    return (u);
}

void random_step(ParamsModel *pm, double step)
{
    unsigned i = 0;
    // calculate the transformed step for all parameters
    // the step is between 0 and max - min
    // there is no reason to transform the original step around one
    // for the parameters that have the transformation around one
    // because the step and the parameters do not have the same meaning

    if (pm->r_flag == TRUE)
    {
        // r[0] should always be 1
        for (i = 1; i < pm->no_groups; i++)
        {
            pm->r[i] = one_step(
                            untran_step(step, pm->r_max - pm->r_min, pm->k),
                            pm->r[i], pm->r_min, pm->r_max);
        }
    }
    if (pm->lambda_flag == TRUE)
    {
        pm->lambda = one_step(
                        untran_step(step, pm->lambda_max - pm->lambda_min, pm->k),
                        pm->lambda, pm->lambda_min, pm->lambda_max);
    }
    if (pm->theta_bar_flag == TRUE)
    {
        pm->theta_bar = one_step(
                        untran_step(step, pm->theta_bar_max - pm->theta_bar_min, pm->k),
                        pm->theta_bar, pm->theta_bar_min, pm->theta_bar_max);
    }
    if (pm->a_flag == TRUE)
    {
        pm->a = one_step(
                        untran_step(step, pm->a_max - pm->a_min, pm->k),
                        pm->a, pm->a_min, pm->a_max);
    }
    if (pm->eps_an_flag == TRUE)
    {
        pm->eps_an = one_step(
                        untran_step(step, pm->eps_an_max - pm->eps_an_min, pm->k),
                        pm->eps_an, pm->eps_an_min, pm->eps_an_max);
    }
    if (pm->eps_cont_flag == TRUE)
    {
        pm->eps_cont = one_step(
                        untran_step(step, pm->eps_cont_max - pm->eps_cont_min, pm->k),
                        pm->eps_cont, pm->eps_cont_min, pm->eps_cont_max);
    }

    if (pm->model < 4)
    {
        for (i = 0; i < pm->no_sel; i++)
        {
            if (pm->sel_flag[i] == FALSE)
            {
                // skip this
                continue;
            }
            pm->sel_params[i] = one_step(
                            untran_step(step,  pm->sel_max[i] - pm->sel_min[i], pm->k),
                            pm->sel_params[i], pm->sel_min[i], pm->sel_max[i]);
        }
    }
    else
    {
        // for model D I have a special randomization of parameters
        // due to the constraint that probabilities need to sum to 1
        double sum = 1;
        double aux_step = 0;
        double max = 0;
        for (i = 0; i < pm->no_sel/2; i++)
        {
            // only randomize if this parameter is estimated
            if (pm->sel_flag[2 * i + 1] == TRUE)
            {
                max = sum < pm->sel_max[2 * i + 1] ?
                                sum : pm->sel_max[2 * i + 1];
                // divide max by the remaining number of parameters
                max /= (pm->no_sel / 2 - (i - 1));
                aux_step = untran_step(step, max - pm->sel_min[2 * i + 1],
                                       pm->k);
                // need to check that pm->sel_params[2*i+1] < max
                pm->sel_params[2 * i + 1] =
                                pm->sel_params[2 * i + 1] < max ?
                                                pm->sel_params[2 * i + 1] : max;
                pm->sel_params[2 * i + 1] = one_step(aux_step,
                                                     pm->sel_params[2 * i + 1],
                                                     pm->sel_min[2 * i + 1],
                                                     max);
            }
            // update the remaining value
            if (pm->sel_flag[2 * i + 1] != 2)
            {
                sum -= pm->sel_params[2*i+1];
            }
        }
        pm->sel_params[pm->sel_fixed] = sum;
    }
}
