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
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#include "basinhopping.h"
#include "likelihood.h"
#include "localoptim.h"
#include "parse.h"
#include "simulate.h"

#define EXIT_HELP 100

int handle_status(int status, char *filename, int param)
{
    switch (status)
    {
        case EXIT_SUCCESS:
        {
            return (EXIT_SUCCESS);
        }
        case FILE_NOT_FOUND:
        {
            fprintf(stderr, "File %s was not found.\n", filename);
            return (EXIT_FAILURE);
        }
        case PARAM_NOT_FOUND:
        {
            fprintf(stderr, "Entry %d was not found in file %s.\n", param,
                    filename);
            return (EXIT_FAILURE);
        }
        default:
        {
            if (status < WRONG_RANGE)
            {
                if (status <= LINE_SHORT)
                {
                    fprintf(stderr,
                            "Line %d in file %s contains too few numbers.\n",
                            LINE_SHORT - status,
                            filename);
                }
                if (status >= LINE_LONG)
                {
                    fprintf(stderr,
                            "Line %d in file %s contains too many numbers.\n",
                            status - LINE_LONG, filename);
                }
            }
            else
            {
                // found wrong range
                // don't do anything
                // as the ranges are updated automatically
//                int wrong = status / WRONG_RANGE;
//                int line_no = status - WRONG_RANGE * wrong;
//                fprintf(stderr,
//                        "Line %d in file %s contains %d numbers that are "
//                        "either invalid or outside permitted limits.\n",
//                        line_no, filename, wrong);
            }
            return (EXIT_FAILURE);
        }
    }

    return (EXIT_SUCCESS);
}

void print_usage(char *argv, FILE *f)
{
    fprintf(f, "%s v1.1\nUsage: %s -d data_file [-m model(A, B, C, D) [K]] [-t]\n\t"
            "{-s m_neut L_neut m_sel L_sel n -i init_file ID"
            " || \n\t[-o optim(bfgs, conj_pr, conj_fr)] [-k kind(s, j, s+j)] [-r range_file ID]"
            "\n\t\t[-i init_file ID [-j]] [-e] [-w] [-g grouping_file ID]"
            "\n\t\t[-p optim_file ID] [-b [basinhop_file ID]]"
            "\n\t\t[-l min] [-v verbose(0, 1, frequency)]}\n",
            argv, argv);
}

int parse_options(int argc, char **argv, unsigned *model, unsigned int *kind,
                  unsigned *method, int *time, int *verbose, int* minutes,
                  int *initial_estimation, int *initial_values,
                  int *div, char **data_file,
                  char **group_file, int *group_id, char **optim_file,
                  int *optim_id, char **range_file, int *range_id,
                  char **init_file, int *init_id, char **basinhop_file,
                  int *basinhop_id, char **cubature_file, int *cubature_id,
                  double **sim)
{
    char opt;
    unsigned i;

    // flags for each required options
    int no_req_flags = 4;
    int req_flags[6] = { 0,   0,   -1,  -1,   -1,   0 };
    char req_opts[6] = { 'A', 'd', 'm', 'o', 'r', 'i' };

    while ((opt = getopt(argc, argv, "d:m:k:o:r:i:p:bs:g:tv:ejwhl:")) != -1)
    {
        switch (opt)
        {
            case 'd':
            {
                (*data_file) = optarg;
                req_flags[1] = 1;
                break;
            }
            case 'm':
            {
                if (strcmp(optarg, "A") == 0)
                {
                    *model = 1;
                }
                else if (strcmp(optarg, "B") == 0)
                {
                    *model = 2;
                }
                else if (strcmp(optarg, "C") == 0)
                {
                    *model = 3;
                }
                else if (strcmp(optarg, "D") == 0)
                {
                    *model = 4;
                    if (optind < argc && *argv[optind] != '-')
                    {
                        *model += atoi(argv[optind]);
                        optind++;
                    }
                    else
                    {
                        fprintf(stderr, "Option -m requires number of selection "
                                "classes K when using model D\n\n");
                        req_flags[0] = 1;
                    }
                }
                else
                {
                    fprintf(stderr, "Invalid argument '%s': -m option "
                            "requires one of these arguments: A, B, C\n\n",
                            optarg);
                    req_flags[0] = 1;
                }
                req_flags[2] = 1;
                break;
            }
            case 'k':
            {
                if (strcmp(optarg, "s") == 0)
                {
                    *kind = 0;
                }
                else if (strcmp(optarg, "s+j") == 0)
                {
                    *kind = 1;
                }
                else if (strcmp(optarg, "j") == 0)
                {
                    *kind = 2;
                }
                else
                {
                    *kind = 2;
                    fprintf(stderr, "Invalid argument '%s': -k option "
                            "requires one of these arguments: s, s+j, j.\n\n",
                            optarg);
                    req_flags[0] = 1;
                }
                break;
            }
            case 'i':
            {
                (*init_file) = optarg;
                if (optind < argc && *argv[optind] != '-')
                {
                    (*init_id) = atoi(argv[optind]);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -i requires TWO arguments "
                            "init_file ID\n\n");
                    req_flags[0] = 1;
                }
                req_flags[5] = 1;
                break;
            }
            case 'o':
            {
                if (strcmp(optarg, "bfgs") == 0)
                {
                    *method = 1;
                }
                else if (strcmp(optarg, "conj_pr") == 0)
                {
                    *method = 2;
                }
                else if (strcmp(optarg, "conj_fr") == 0)
                {
                    *method = 3;
                }
                else
                {
                    fprintf(stderr, "Invalid argument '%s': -o' option "
                            "requires one of these arguments: bfgs, "
                            "conj_pr, conj_fr \n\n",
                            optarg);
                    req_flags[0] = 1;
                }
                req_flags[3] = 1;

                break;
            }
            case 'r':
            {
                (*range_file) = optarg;
                if (optind < argc && *argv[optind] != '-')
                {
                    (*range_id) = atoi(argv[optind]);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -r requires TWO arguments "
                            "range_file ID\n\n");
                    req_flags[0] = 1;
                }
                req_flags[4] = 1;
                break;
            }
            case 'p':
            {
                (*optim_file) = optarg;
                if (optind < argc && *argv[optind] != '-')
                {
                    (*optim_id) = atoi(argv[optind]);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -p requires TWO arguments "
                            "param_file ID\n\n");
                    req_flags[0] = 1;
                }
                break;
            }
            case 'b':
            {
                (*basinhop_id) = -2;
                if (optind < argc && *argv[optind] != '-')
                {
                    (*basinhop_file) = argv[optind];
                    optind++;
                    if (optind < argc && *argv[optind] != '-')
                    {
                        (*basinhop_id) = atoi(argv[optind]);
                        optind++;
                    }
                    else
                    {
                        fprintf(stderr, "Option -b requires TWO arguments "
                                "basin_file ID\n\n");
                        req_flags[0] = 1;
                    }
                }
                break;
            }
            case 's':
            {
                (*sim)[0] = atoi(optarg);
                i = 1;
                while (optind < argc && *argv[optind] != '-' && i < 5)
                {
                    (*sim)[i] = atoi(argv[optind]);
                    optind++;
                    i++;
                }
                if (i != 5)
                {
                    fprintf(stderr, "Option -s requires FIVE arguments "
                            "m_neut L_neut m_sel L_sel n\n\n");
                    req_flags[0] = 1;
                }
                if (req_flags[3] != -1)
                {
                    req_flags[3] += 2;
                } else
                {
                    req_flags[3] = 2;
                }

                break;
            }
            case 'g':
            {
                (*group_file) = optarg;
                if (optind < argc && *argv[optind] != '-')
                {
                    (*group_id) = atoi(argv[optind]);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -g requires TWO arguments "
                            "grouping_file ID\n\n");
                    req_flags[0] = 1;
                }
                break;
            }
            case 't':
            {
                *time = TRUE;
                break;
            }
            case 'v':
            {
                *verbose = atoi(optarg);
                break;
            }
            case 'l':
            {
                *minutes = atoi(optarg);
                break;
            }
            case 'e':
            {
                *initial_estimation = TRUE;
                break;
            }
            case 'j':
            {
                *initial_estimation = TRUE;
                *initial_values = TRUE;
                break;
            }
            case 'w':
            {
                *div = FALSE;
                break;
            }
            case 'h':
            {
                print_usage(argv[0], stdout);
                return (EXIT_HELP);
                break;
            }
            default:
            {
                print_usage(argv[0], stderr);
                return (EXIT_FAILURE);
            }
        }
    }

//    // if -m is given and no -s, I want -r
//    if (req_flags[2] != -1 && req_flags[3] < 2)
//    {
//        if (req_flags[4] == -1)
//        {
//            fprintf(stderr, "Option -m requires option -r\n\n");
//            req_flags[0] = 1;
//        }
//    }

    // verify and generate any errors reading the arguments
    for (i = 1; i < no_req_flags - 1; i++)
    {
        if (req_flags[i] == 0)
        {
            fprintf(stderr, "Option -%c is missing\n\n", req_opts[i]);
            req_flags[0] = 1;
        }
    }

    if (req_flags[3] == 1)
    {
        if (req_flags[4] == 0)
        {
            fprintf(stderr, "Option -r is missing "
                    "and it is required by -o\n\n");
            req_flags[0] = 1;
        }
    }
    else if (req_flags[3] == 2)
    {
        if (req_flags[5] == 0)
        {
            fprintf(stderr, "Option -i is missing "
                    "and it is required by -s\n\n");
            req_flags[0] = 1;
        }
    }
    else if (req_flags[3] != -1)
    {
        fprintf(stderr, "Options -o and -s cannot be used together\n\n");
        req_flags[0] = 1;
    }

    // verify if init file is given for discretized DFE
    if (*model > 3 && (*init_id) == -1)
    {
        fprintf(stderr,"Option -i is missing "
                "and it is required by -m D\n\n");
        req_flags[0] = 1;
    }

    // verify if init file is given when using -j
    if (*initial_values == TRUE && (*init_id) == -1)
    {
        fprintf(stderr, "Option -i is missing "
                "and it is required by -j\n\n");
        req_flags[0] = 1;
    }

    if (req_flags[0] == 1)
    {
        print_usage(argv[0], stderr);
        return (EXIT_FAILURE);
    }

    if (*init_id == -1)
    {
        // use automatic initial estimation
        *initial_estimation = TRUE;
    }

    return (EXIT_SUCCESS);
}

void set_model(ParamsModel *pm, int model)
{
    if (model < 4)
    {
        pm->model = model;
    }
    else
    {
        pm->model = 4;
    }

    switch (pm->model)
    {
        case 1:
        {
            pm->no_sel = 3;
            break;
        }
        case 2:
        {
            pm->no_sel = 4;
            break;
        }
        case 3:
        {
            pm->no_sel = 4;
            break;
        }
        case 4:
        {
            // I need twice as many
            // for both probabilities and selection coefficients
            // probabilities will be in odd positions
            pm->no_sel = 2 * (model - 4);
            break;
        }
    }
}

int simulate(unsigned model, double *sim, char *data_file, char *init_file,
             int init_id, char *cubature_file, int cubature_id)
{
    int parsed = EXIT_SUCCESS;
    int status = EXIT_SUCCESS;
    char c;

    fprintf(stderr, "Warning: simulating data, %s will be overwritten.\n",
            data_file);
    fflush(stderr);
    printf("Do you want to continue with simulation? (Y/N): ");
    fflush(stdout);
    c = getchar();
    if (c != 'Y' && c != 'N')
    {
        printf("Unknown input %c. Terminating...\n\n", c);
        return (EXIT_SUCCESS);
    }
    if (c != 'Y')
    {
        printf("You chose no. Terminating...\n\n");
        return (EXIT_SUCCESS);
    }

    // initialize the parameters
    ParamsModel pm;
    initialize_params_model(&pm);
    set_model(&pm, model);
    pm.n = sim[4];
    allocate_selection_params(&pm);
    allocate_grouping(&pm, NULL);

    parsed = parse(5, init_file, init_id, &pm);
    status += handle_status(parsed, init_file, init_id);

    if (status != EXIT_SUCCESS)
    {
        free_params_model(&pm);
        return (EXIT_FAILURE);
    }

    status = sim_data(&pm, sim, data_file);
    free_params_model(&pm);

    if (status != EXIT_SUCCESS)
    {
        return (EXIT_FAILURE);
    }

    return (EXIT_SUCCESS);
}

int optimize(unsigned model, unsigned kind, unsigned method, char *data_file,
             char *group_file, int group_id, char *optim_file, int optim_id,
             char *range_file, int range_id, char *init_file, int init_id,
             char *basinhop_file, int basinhop_id, char *cubature_file,
             int cubature_id, int verbose, int minutes,  int initial_estimation,
             int initial_values, int div)
{
    int status = EXIT_SUCCESS;
    int parsed = EXIT_SUCCESS;

    Params p;
    initialize_params(&p);
    set_model(p.pm, model);
    // set divergence info
    p.pm->div_flag = div;
    if (p.pm->div_flag == FALSE)
    {
        // do not estimate lambda
        p.pm->lambda_flag = FALSE;
    }
    p.kind = kind;
    p.verbose = verbose;

    ParamsOptim po;
    initialize_params_optim(&po);
    po.minutes = minutes;

    ParamsBasinHop pb;
    initialize_params_basin_hop(&pb);
    p.pm->inital_estimation = initial_estimation;
    p.pm->initial_values = initial_values;

    // the params object is allocated in parse_data
    parsed = parse_data(data_file, &p, &p.pm->n);
    status += handle_status(parsed, data_file, 0);
    if (status != EXIT_SUCCESS)
    {
        free_params(&p);
        return (EXIT_FAILURE);
    }
    printf("---- Performing inference on %s using model ", data_file);
    switch (model)
    {
        case 1:
        {
            printf("A\n\n");
            break;
        }
        case 2:
        {
            printf("B\n\n");
            break;
        }
        case 3:
        {
            printf("C\n\n");
            break;
        }
        case 4:
        {
            printf("D\n\n");
            break;
        }
    }

    parsed = parse(6, group_file, group_id, p.pm);
    status += handle_status(parsed, group_file, group_id);

    if (range_file != NULL)
    {
        parsed = parse(0, range_file, range_id, p.pm);
        status += handle_status(parsed, range_file, range_id);
    }

    parsed = parse(2, optim_file, optim_id, &po);
    status += handle_status(parsed, optim_file, optim_id);

    parsed = parse(3, basinhop_file, basinhop_id, &pb);
    status += handle_status(parsed, basinhop_file, basinhop_id);

    initialize_selection_params(p.pm, range_file);
    if (init_id != -1)
    {
        // if I don't use initial_estimation, make sure initial_values is TRUE
        if (p.pm->inital_estimation == FALSE)
        {
            p.pm->initial_values = TRUE;
        }
        parsed = parse(1, init_file, init_id, p.pm);
        status += handle_status(parsed, init_file, init_id);
    }


    if (status != EXIT_SUCCESS)
    {
        free_params(&p);
        return (EXIT_FAILURE);
    }

    allocate_params_basin_hop(&pb, p.pm);

    // define optimization algorithm
    const gsl_multimin_fdfminimizer_type *type = NULL;
    gsl_set_error_handler_off();

    switch (method)
    {
        case 1:
        {
            type = gsl_multimin_fdfminimizer_vector_bfgs2;
            break;
        }
        case 2:
        {
            type = gsl_multimin_fdfminimizer_conjugate_pr;
            break;
        }
        case 3:
        {
            type = gsl_multimin_fdfminimizer_conjugate_fr;
            break;
        }
    }

    if (p.counts_neut[0].sfs[p.pm->n - 1] == -1)
    {
        printf("---- Data does not contain divergence counts. Using option -w.\n");
        // make sure that div and lambda are not used
        p.pm->div_flag = FALSE;
        p.pm->lambda_flag = FALSE;
    }

    status = run_basin_hopping(&pb, type, po, &p);
    // free memory
    free_params(&p);
    free_params_basin_hop(&pb);

    if (status == WRONG_INIT)
    {
        return (EXIT_FAILURE);
    }
    return (EXIT_SUCCESS);
}

int main(int argc, char **argv)
{
    // this is for Eclipse working with stdout and stderr
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    char *data_file = NULL;
    char *group_file = NULL;
    char *optim_file = NULL;
    char *range_file = NULL;
    char *init_file = NULL;
    char *basinhop_file = NULL;
    char *cubature_file = NULL;

    int group_id = -1;
    int optim_id = -1;
    int range_id = -1;
    int init_id = -1;
    int basinhop_id = -1;
    int cubature_id = -1;
    int initial_estimation = FALSE;
    int initial_values = FALSE;
    int div = TRUE; // use divergence data

    double *sim = malloc(sizeof(double) * 5);
    sim[0] = -1;

    // default model: 3 (C)
    unsigned model = 3;
    // default kind: 2 (joint likelihood)
    unsigned kind = 2;
    // default method: 1 (BFGS)
    unsigned method = 1;
    int time = FALSE;
    int verbose = 0;
    // default limit on running time should be 5h - 300 minutes!
    int minutes = 300;

    int status = parse_options(argc, argv, &model, &kind, &method, &time, &verbose,
                               &minutes, &initial_estimation, &initial_values,
                               &div, &data_file,
                               &group_file, &group_id, &optim_file, &optim_id,
                               &range_file, &range_id, &init_file, &init_id,
                               &basinhop_file, &basinhop_id, &cubature_file,
                               &cubature_id, &sim);

    if (status != EXIT_SUCCESS)
    {
        if (sim)
        {
            free(sim);
        }

        if (status == EXIT_FAILURE)
        {
            printf("\n");
            return (EXIT_FAILURE);
        }
        else
        {
            printf("\n");
            return(EXIT_SUCCESS);
        }
    }

    int i;
    printf("---- Running command\n---- ./polyDFE ");
    for (i = 1; i < argc; i++)
    {
        printf("%s ",argv[i]);
    }
    printf("\n");

    // turn off gsl error handler
    gsl_set_error_handler_off();

    clock_t begin = 0, end;
    double time_spent;
    if (time == TRUE)
    {
        begin = clock();
    }

    if (sim[0] != -1)
    {
        status = simulate(model, sim, data_file, init_file, init_id,
                          cubature_file, cubature_id);
    }
    else
    {
        status = optimize(model, kind, method, data_file, group_file, group_id,
                          optim_file, optim_id, range_file, range_id, init_file,
                          init_id, basinhop_file, basinhop_id, cubature_file,
                          cubature_id, verbose, minutes,
                          initial_estimation, initial_values, div);
    }

    if (time == TRUE && status != EXIT_FAILURE)
    {
        end = clock();
        time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
        int minutes = floor(time_spent / 60);
        int seconds = floor(time_spent - minutes * 60);
        printf("---- Running time: %d minutes and %d seconds\n\n\n", minutes,
               seconds);
    }

    free(sim);

    return (status);
}
