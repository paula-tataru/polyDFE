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

#include "parse.h"

#include <ctype.h>
#include <float.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "basinhopping.h"
#include "localoptim.h"

// getline function for systems that don't have it in the default libraries
#if !defined(_SYS_STDIO_H_) && !defined(__USE_XOPEN2K8)
ssize_t getline(char **lineptr, size_t *n, FILE *stream)
{
    char *bufptr = NULL;
    char *p = bufptr;
    size_t size;
    int c;
    if (lineptr == NULL)
    {
        return (-1);
    }
    if (stream == NULL)
    {
        return (-1);
    }
    if (n == NULL)
    {
        return (-1);
    }
    bufptr = *lineptr;
    size = *n;
    c = fgetc(stream);
    if (c == EOF)
    {
        return (-1);
    }
    if (bufptr == NULL)
    {
        bufptr = malloc(128);
        if (bufptr == NULL)
        {
            return (-1);
        }
        size = 128;
    }
    p = bufptr;
    while (c != EOF)
    {
        if ((p - bufptr) > (size - 1))
        {
            size = size + 128;
            bufptr = realloc(bufptr, size);
            if (bufptr == NULL)
            {
                return (-1);
            }
        }
        *p++ = c;
        if (c == '\n')
        {
            break;
        }
        c = fgetc(stream);
    }
    *p++ = '\0';
    *lineptr = bufptr;
    *n = size;
    return (p - bufptr - 1);
}

#endif

char *trim(char *str)
{
    // trim leading space
    while (isspace((int ) *str))
    {
        str++;
    }

    return (str);
}

int parse_line(char *buf, unsigned param_len, double *temp)
{
    unsigned i;
    char *left;

    for (i = 0; i < param_len; i++)
    {
        // strip for spaces
        left = trim(buf);
        if (left[0] == '\0')
        {
            return (LINE_SHORT);
        }
        temp[i] = strtod(buf, &buf);
    }

    // strip for spaces
    left = trim(buf);
    if (left[0] != '\0')
    {
        return (LINE_LONG);
    }

    return (EXIT_SUCCESS);
}

int handle_line_status(int status, int line)
{
    switch (status)
    {
        case EXIT_SUCCESS:
        {
            return (EXIT_SUCCESS);
        }
        case LINE_SHORT:
        {
            return (LINE_SHORT - line);
        }
        case LINE_LONG:
        {
            return (LINE_LONG + line);
        }
        default:
        {
            // found wrong range
            return (status + line);
        }
    }

    return (EXIT_SUCCESS);
}

int parse_data(char *filename, Params *p, unsigned *n)
{
    FILE *file_in = fopen(filename, "r");
    if (file_in == NULL)
    {
        return (FILE_NOT_FOUND);
    }

    char *buf = NULL;
    size_t buf_len = 0;
    unsigned i = -1;
    unsigned j = 0;
    unsigned line_cnt = 0;
    int pp;
    int status = EXIT_SUCCESS;
    int no_fields = 0;
    double *temp = NULL;

    while ((pp = getline(&buf, &buf_len, file_in)) != -1)
    {
        line_cnt++;
        // ignore empty line and line beginning with '#'
        if (pp <= 1 || buf[0] == '#' || buf[0] == '\n'|| buf[0] == '\r')
        {
            continue;
        }

        if (i == -1)
        {
            // first line contains no of fragments and n
            int input_len = 3;
            temp = malloc(sizeof(double) * input_len);
            status = parse_line(buf, input_len, temp);
            if (status != EXIT_SUCCESS)
            {
                break;
            }

            // verify the validity of len and n, stored in temp
            // they should be at least 1
            status = EXIT_SUCCESS;
            for (j = 0; j < input_len; j++)
            {
                status += check_lim(temp[j], 1, DBL_MAX);
            }
            if (status != EXIT_SUCCESS)
            {
                break;
            }

            p->no_neut = (unsigned) temp[0];
            p->no_sel = (unsigned) temp[1];
            (*n) = (unsigned) temp[2];

            allocate_params(p);
            // n - 1 SFS + 1 divergence + 2 lengths
            no_fields = (*n) + 2;

            free(temp);
            temp = malloc(sizeof(double) * no_fields);
        }
        else
        {
            status = parse_line(buf, no_fields, temp);
            if (status != EXIT_SUCCESS)
            {
                // check if maybe divergence data is missing
            	if (status == LINE_SHORT)
            	{
            		// n - 1 SFS + 1 length
            		no_fields = (*n);
            		// re-read the line to see if it matches
            		status = parse_line(buf, no_fields, temp);
            		if (status == EXIT_SUCCESS)
            		{
            			// mark no divergence in the data
            			p->counts_neut[0].sfs[(*n) - 1] = -1;
            			p->pm->div_flag = FALSE;
            			p->pm->lambda_flag = FALSE;
            		}
            		else
            		{
            			break;
            		}
            	}
            }

            // verify the validity of counts stored in temp
            // they should be at least 0
            for (j = 0; j < no_fields; j++)
            {
                status += check_lim(temp[j], 0, DBL_MAX);
            }
            if (status != EXIT_SUCCESS)
            {
                break;
            }

            if (i < p->no_neut)
            {
                // first sfs
                for (j = 0; j < (*n) - 1; j++)
                {
                    p->counts_neut[i].sfs[j] = temp[j];
                }
                j = (*n) - 1;
                p->counts_neut[i].len_poly = temp[j];

                if (no_fields > (*n))
                {
                    // now divergence
                    p->counts_neut[i].sfs[j] = temp[j + 1];
                    p->counts_neut[i].len_div = temp[j + 2];
                }
            }
            else
            {
                // first sfs
                for (j = 0; j < (*n) - 1; j++)
                {
                    p->counts_sel[i - p->no_neut].sfs[j] = temp[j];
                }
                j = (*n) - 1;
                p->counts_sel[i - p->no_neut].len_poly = temp[j];

                if (no_fields > (*n))
                {
                    // now divergence
                    p->counts_sel[i - p->no_neut].sfs[j] = temp[j + 1];
                    p->counts_sel[i - p->no_neut].len_div = temp[j + 2];
                }
            }
        }

        i++;
    }

    fclose(file_in);
    free(buf);
    free(temp);

    if (status != EXIT_SUCCESS)
    {
        return (handle_line_status(status, line_cnt));
    }

    return (EXIT_SUCCESS);
}

int parse_range(ParamsModel *pm, double *temp)
{
    int status = EXIT_SUCCESS;
    unsigned i = 1;
    unsigned j = 0;

    // I do not use eps_cont anymore
    // but I still read it for backward compatibility
    double eps_cont = 0;

    pm->k = temp[i++];
    pm->eps_an_min = temp[i++];
    pm->eps_an_max = temp[i++];
    eps_cont = temp[i++];
    eps_cont = temp[i++];
    pm->lambda_min = temp[i++];
    pm->lambda_max = temp[i++];

    // get rid of warning for eps_cont
    eps_cont++;

    // parameters for mutation distribution
    pm->theta_bar_min = temp[i++];
    pm->theta_bar_max = temp[i++];
    pm->a_min = temp[i++];
    pm->a_max = temp[i++];

    // parameters for DFE
    if (pm->model < 4)
    {
        for (j = 0; j < pm->no_sel; j++)
        {
            pm->sel_min[j] = temp[i++];
            pm->sel_max[j] = temp[i++];
        }
    }
    else
    {
        // if using the discretized DFE, I have just one range for everything
        pm->sel_min[0] = temp[i++];
        pm->sel_max[0] = temp[i++];
        for (j = 1; j < pm->no_sel; j++)
        {
            pm->sel_min[j] = pm->sel_min[0];
            pm->sel_max[j] = pm->sel_max[0];
        }
    }

    pm->r_min = temp[i++];
    pm->r_max = temp[i++];

    // verify that I have sensible input
    status += check_lim(pm->k, 0, DBL_MAX);
    status += check_lim(pm->eps_an_min, 0, pm->eps_an_max);
    status += check_lim(pm->eps_an_max, 0, 1);
    status += check_lim(pm->lambda_min, 0, pm->lambda_max);
    status += check_lim(pm->lambda_max, 0, 10);
    status += check_lim(pm->theta_bar_min, 0, pm->theta_bar_max);
    status += check_lim(pm->theta_bar_max, 0, DBL_MAX);
    status += check_lim(pm->a_min, 0, pm->a_max);
    status += check_lim(pm->a_max, 0, DBL_MAX);

    // parameters for DFE
    // first, check that the minimum is smaller than the maximum
    // and that the values are found between -DBL_MAX and DBL_MAX
    for (i = 0; i < pm->no_sel; i++)
    {
        status += check_lim(pm->sel_min[i], -DBL_MAX, pm->sel_max[i]);
        status += check_lim(pm->sel_max[i], pm->sel_min[i], DBL_MAX);
    }

    if (pm->model != 4)
    {
        // s_bar or s_d
        status += check_lim(pm->sel_max[0], -DBL_MAX, 0);
        // b
        status += check_lim(pm->sel_min[1], 0, pm->sel_max[1]);
        if (pm->model == 1)
        {
            // s_max
            status += check_lim(pm->sel_min[2], 0, pm->sel_max[2]);
        }
        else
        {
            // p_b
            status += check_lim(pm->sel_min[2], 0, pm->sel_max[2]);
            status += check_lim(pm->sel_max[2], 0, 1);
            // s_b
            status += check_lim(pm->sel_min[3], 0, pm->sel_max[3]);
        }
    }
    else
    {
        // all the limits are the same - enough to check for entry 0
        status += check_lim(pm->sel_min[0], 0, pm->sel_max[0]);
        status += check_lim(pm->sel_max[0], 0, 1);
    }

    status += check_lim(pm->r_min, 0, pm->r_max);
    status += check_lim(pm->r_max, 0, DBL_MAX);

    return (status);
}

int parse_init(ParamsModel *pm, double *temp, int check_range, int initial_values)
{
    int status = EXIT_SUCCESS;
    unsigned i = 1;
    unsigned j = 0;

    // I do not use eps_cont anymore
    // but I still read it for backward compatibility
    double eps_cont = 0;

    // 1 is TRUE and 0 is FALSE but
    // 1 is fix and 0 is estimate
    // so I need to get 1 - temp!
    // flag to share a parameter will then be -1
    // as I set it to 2 in the init file
    pm->eps_an_flag = 1 - temp[i++];
    pm->eps_an = temp[i++];
    eps_cont = 1 - temp[i++];
    eps_cont = temp[i++];
    pm->lambda_flag = 1 - temp[i++];
    pm->lambda = temp[i++];

    // get rid of warning for eps_cont
    eps_cont++;

    // parameters for mutation distribution
    pm->theta_bar_flag = 1 - temp[i++];
    pm->theta_bar = temp[i++];
    pm->a_flag = 1 - temp[i++];
    pm->a = temp[i++];

    // parameters for DFE
    if (pm->model < 4)
    {
        for (j = 0; j < pm->no_sel; j++)
        {
            pm->sel_flag[j] = 1 - temp[i++];
            pm->sel_params[j] = temp[i++];
        }
    }
    else
    {
        // set the flag of the first parameter that is estimated
        // to be FALSE, as this parameter should be 1 - sum(rest)
        int found_true = FALSE;
        for (j = 0; j < pm->no_sel/2; j++)
        {
            // set the flag of the coefficients to FALSE!
            pm->sel_flag[2*j] = FALSE;
            pm->sel_flag[2*j+1] = 1 - temp[i++];
            pm->sel_params[2*j] = temp[i++];
            pm->sel_params[2*j+1] = temp[i++];

            if (pm->sel_flag[2*j+1] == TRUE && found_true == FALSE)
            {
                found_true = TRUE;
                // mark it by setting it to something else than FALSE
                pm->sel_flag[2*j+1] = 2;
            }
        }
    }

    // r[0] should always be 1
    pm->r_flag = 1 - temp[i++];
    pm->r[0] = 1;
    for (j = 1; j < pm->no_groups; j++)
    {
        pm->r[j] = temp[i++];
    }

    if (check_range == TRUE)
    {
        // verify that parameters are in the correct range
        // I only want to verify the provided init values is within the range
        // if the parameter is used as it is
        // (not estimated because of the flag or used as init)
        status += check_lim(pm->eps_an_flag, -1, 1);
        if (pm->eps_an_flag == FALSE || initial_values == TRUE)
        {
            status += check_lim_update(&pm->eps_an, &pm->eps_an_min, &pm->eps_an_max,
                                               pm->eps_an_flag, "eps_an");
        }
        status += check_lim(pm->lambda_flag, -1, 1);
        if (pm->lambda_flag == FALSE || initial_values == TRUE)
        {
            status += check_lim_update(&pm->lambda, &pm->lambda_min,
                                       &pm->lambda_max, pm->lambda_flag,
                                       "lambda");
        }
        status += check_lim(pm->theta_bar_flag, -1, 1);
        if (pm->theta_bar_flag == FALSE || initial_values == TRUE)
        {
            status += check_lim_update(&pm->theta_bar, &pm->theta_bar_min,
                                       &pm->theta_bar_max, pm->theta_bar_flag,
                                       "theta_bar");
        }
        status += check_lim(pm->a_flag, -1, 1);
        // allow for a to be -1 and outside range
        if (pm->a != -1)
        {
            if (pm->a_flag == FALSE || initial_values == TRUE)
            {
                status += check_lim_update(&pm->a, &pm->a_min, &pm->a_max,
                                           pm->a_flag, "a");
            }
        }

        // parameters for DFE
        if (pm->model < 4)
        {
            for (i = 0; i < pm->no_sel; i++)
            {
                status += check_lim(pm->sel_flag[i], -1, 1);
                if (pm->sel_flag[i] == FALSE || initial_values == TRUE)
                {
                    status += check_lim_update(&pm->sel_params[i],
                                               &pm->sel_min[i], &pm->sel_max[i],
                                               pm->sel_flag[i],
                                               "a selection parameter");
                }
            }
        }
        else
        {
            // for discretized DFE, check only the probabilities
            // and also make sure they sum to 1
            double sum = 0;
            for (i = 0; i < pm->no_sel/2; i++)
            {
                if (2*i+1 != pm->sel_fixed)
                {
                    status += check_lim(pm->sel_flag[2 * i + 1], 0, 1);
                }
                if (pm->sel_flag[2 * i + 1] == FALSE
                                || initial_values == TRUE)
                {
                    status += check_lim_update(&pm->sel_params[2 * i + 1],
                                               &pm->sel_min[2 * i + 1],
                                               &pm->sel_max[2 * i + 1],
                                               pm->sel_flag[2 * i + 1],
                                               "a selection parameter");
                }
                sum += pm->sel_params[2 * i + 1];
            }
            status += check_lim(sum, 1, 1);
        }

        status += check_lim(pm->r_flag, -1, 1);
        for (j = 1; j < pm->no_groups; j++)
        {
            if (pm->r_flag == FALSE || initial_values == TRUE)
            {
                status += check_lim_update(&pm->r[j], &pm->r_min, &pm->r_max,
                                           pm->r_flag, "r");
            }
        }
    }

    return (status);
}

int parse_optim(ParamsOptim *po, double *temp)
{
    int status = EXIT_SUCCESS;
    unsigned i = 1;

    po->eps_abs = temp[i++];
    po->step_size = temp[i++];
    po->tol = temp[i++];
    po->max_iter = temp[i++];

    // verify that I have sensible input
    status += check_lim(po->eps_abs, 0, 1);
    status += check_lim(po->step_size, 0, DBL_MAX);
    status += check_lim(po->tol, 0, 1);
    status += check_lim(po->max_iter, 1, DBL_MAX);

    return (status);
}

int parse_basinhop(ParamsBasinHop *pb, double *temp)
{
    int status = EXIT_SUCCESS;
    unsigned i = 1;

    pb->max_same = temp[i++];
    pb->max_iter = temp[i++];
    pb->temp = temp[i++];
    pb->step = temp[i++];
    pb->accept_rate = temp[i++];
    pb->interval = temp[i++];
    pb->factor = temp[i++];

    // verify that I have sensible input
    status += check_lim(pb->max_same, 1, DBL_MAX);
    status += check_lim(pb->max_iter, 0, DBL_MAX);
    status += check_lim(pb->temp, 0, DBL_MAX);
    status += check_lim(pb->step, 0, DBL_MAX);
    status += check_lim(pb->accept_rate, 0, 1);
    status += check_lim(pb->interval, 1, DBL_MAX);
    status += check_lim(pb->factor, 0, 1);

    return (status);
}

int no_sel_params(int model, double k)
{
    switch (model)
    {
        case 1:
        {
            return (2 * 3);
            break;
        }
        case 2:
        {
            return (2 * 4);
            break;
        }
        case 3:
        {
            return (2 * 4);
            break;
        }
        case 4:
        {
            return ((int) 3 * k);
            break;
        }
    }
    return (2);
}

int parse(int what, char *filename, int id, void *pv, int initial_values)
{
    if (id == -1 && what != 6)
    {
        // no id and I am not looking for grouping
        // then it's an optional parsing
        return (EXIT_SUCCESS);
    }

    FILE *file_in = NULL;
    if (id >= 0)
    {
        file_in = fopen(filename, "r");
        if (file_in == NULL)
        {
            return (FILE_NOT_FOUND);
        }
    }

    char *buf = NULL;
    size_t buf_len = 0;
    unsigned line_cnt = 0;
    int pp;
    int status = EXIT_SUCCESS;
    int found = -1;
    int no_params = -1;

    switch (what)
    {
        case 0:
        {
            // range
            // non selection part
            no_params = 1 + 1 + 2 * 5 + 2;
            // selection part
            no_params += no_sel_params(((ParamsModel *) pv)->model, 2.0/3);
            break;
        }
        case 1:
        {
            // init
            // non selection part
            no_params = 1 + 2 * 5 + ((ParamsModel *) pv)->no_groups;
            // selection part
            no_params += no_sel_params(((ParamsModel *) pv)->model,
                                            ((ParamsModel *) pv)->no_sel/2);
            break;
        }
        case 2:
        {
            // local optimum
            no_params = 5;
            break;
        }
        case 3:
        {
            // basin hopping
            if (id == -2)
            {
                // default basin hopping
                set_default_max_iter((ParamsBasinHop *) pv);
                printf("---- Running basin hopping with default values\n");
                return (EXIT_SUCCESS);
            }
            no_params = 8;
            break;
        }
        case 4:
        {
            // this used to be for cubature
            break;
        }
        case 5:
        {
            // init for simulation
            // non selection part
            no_params = 1 + 2 * 5 + ((ParamsModel *) pv)->no_groups;
            // selection part
            no_params += no_sel_params(((ParamsModel *) pv)->model,
                                       ((ParamsModel *) pv)->no_sel/2);
            // if no divergence and no grouping, I have one extra r
            if (((ParamsModel *) pv)->div_flag == FALSE
                            && ((ParamsModel *) pv)->no_groups
                                    == ((ParamsModel *) pv)->n - 1)
            {
                no_params++;
            }
            break;
        }
        case 6:
        {
            // grouping
            if (id == -1)
            {
                // no grouping, set default
                allocate_grouping((ParamsModel *) pv, NULL);
                return (EXIT_SUCCESS);
            }
            no_params = ((ParamsModel *) pv)->n;
            break;
        }
    }

    double *temp = malloc(sizeof(double) * no_params);

    while ((pp = getline(&buf, &buf_len, file_in)) != -1)
    {
        line_cnt++;
        // ignore empty lines and lines beginning with '#'
        if (pp <= 1 || buf[0] == '#')
        {
            continue;
        }

        status = parse_line(buf, no_params, temp);
        if (temp[0] == id)
        {
            // if parsing grouping, I have less params than I asked for
            // if parsing init, I might have more params than I asked for
            if (status != EXIT_SUCCESS && !(status == LINE_SHORT && what == 6)
                            && !(status == LINE_LONG && (what == 1 || what == 5)))
            {
                break;
            }
            else
            {
                status = EXIT_SUCCESS;
            }

            found = 0;
            switch (what)
            {
                case 0:
                {
                    // range
                    status = parse_range((ParamsModel *) pv, temp);
                    break;
                }
                case 1:
                {
                    // init
                    status = parse_init((ParamsModel *) pv, temp, TRUE, initial_values);
                    break;
                }
                case 2:
                {
                    // local optimum
                    status = parse_optim((ParamsOptim *) pv, temp);
                    break;
                }
                case 3:
                {
                    // basin hopping
                    status = parse_basinhop((ParamsBasinHop *) pv, temp);
                    break;
                }
                case 4:
                {
                    // this used to be for cubature
                    break;
                }
                case 5:
                {
                    // init for simulation
                    status = parse_init((ParamsModel *) pv, temp, FALSE, initial_values);
                    break;
                }
                case 6:
                {
                    // grouping
                    status = allocate_grouping((ParamsModel *) pv, temp);
                    break;
                }
            }
            break;
        }
        else
        {
            status = EXIT_SUCCESS;
        }
    }

    fclose(file_in);
    free(buf);
    free(temp);

    if (status != EXIT_SUCCESS)
    {
        return (handle_line_status(status, line_cnt));
    }

    if (found == -1)
    {
        return (PARAM_NOT_FOUND);
    }

    return (EXIT_SUCCESS);
}
