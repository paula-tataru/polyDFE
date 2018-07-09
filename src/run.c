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
            fprintf(stderr, "Entry with ID %d was not found in file %s.\n", param,
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
               // don't do anything as the ranges as updated automatically
               int wrong = status / WRONG_RANGE;
               int line_no = status - WRONG_RANGE * wrong;
               fprintf(stderr,
                       "Line %d in file %s contains ", line_no, filename);
               if (wrong == 1)
               {
            	   fprintf(stderr, "1 number that is ");
               }
               else
               {
            	   fprintf(stderr, "%d numbers that are ", wrong);
               }
               fprintf(stderr, "either invalid or outside permitted limits.\n");
            }
            return (EXIT_FAILURE);
        }
    }

    return (EXIT_SUCCESS);
}

void print_usage(char *argv, FILE *f)
{
    fprintf(f, "%sv2.0\nUsage: %s -d data_file_1[:data_file_2:...:data_file_j]"
    		"\n\t[-m model(A, B, C, D) [K]] [-i init_file ID_1[:ID_2:...:ID_j]] [-t]"
            "\n\t{-s m_neut L_neut m_sel L_sel n ||"
            "\n\t[-o optim(bfgs, conj_pr, conj_fr, simplex)] [-k kind(s, j, s+j)]"
            "\n\t\t[-r range_file ID_1[:ID_2:...:ID_j]] [-i init_file ID [-j]]"
            "\n\t\t[-e] [-w] [-g grouping_file ID_1[:ID_2:...:ID_j]] [-p optim_file ID]"
            "\n\t\t[-b [basinhop_file ID]] [-l min] [-v verbose(0, 1, frequency)]}\n",
            argv, argv);
}

// from https://www.quora.com/How-do-you-write-a-C-program-to-split-a-string-by-a-delimiter-I-think-it-uses-ktochar
// updated to use strtok instead of strtok_r
char **strsplit(const char* str, const char* delim, size_t* numtokens) {
    // copy the original string so that we don't overwrite parts of it
    // (don't do this if you don't need to keep the old line,
    // as this is less efficient)
    char *s = strdup(str);
    // these three variables are part of a very common idiom to
    // implement a dynamically-growing array
    size_t tokens_alloc = 1;
    size_t tokens_used = 0;
    char **tokens = calloc(tokens_alloc, sizeof(char*));
    char *token;

    token = strtok(s, delim);
    while (token != NULL)
    {
    	// check if we need to allocate more space for tokens
    	if (tokens_used == tokens_alloc) {
    		tokens_alloc *= 2;
    		tokens = realloc(tokens, tokens_alloc * sizeof(char*));
    	}
    	tokens[tokens_used++] = strdup(token);
    	token = strtok(NULL, delim);
    }
    // cleanup
    if (tokens_used == 0) {
        free(tokens);
        tokens = NULL;
    } else {
        tokens = realloc(tokens, tokens_used * sizeof(char*));
    }
    *numtokens = tokens_used;
    free(s);

    return (tokens);
}

int *strsplit_no(const char* str, const char* delim, size_t* numtokens)
{
    size_t i = 0;
    char **tokens = strsplit(str, delim, numtokens);
    int *tokens_no = malloc(sizeof(int) * *numtokens);

    // need to apply atoi
    // and free original tokens
    for (i = 0; i < *numtokens; i++)
    {
        tokens_no[i] = atoi(tokens[i]);
        free(tokens[i]);
    }
    free(tokens);

    return (tokens_no);
}

int parse_options(int argc, char **argv, unsigned *model, unsigned int *kind,
                  unsigned *method, int *time, int *verbose, int* minutes,
                  int *initial_estimation, int *initial_values,
                  int *div, char ***data_files, size_t *no_files,
                  char **group_file, int **group_id, size_t *no_group_id,
				  char **optim_file, int *optim_id,
				  char **range_file, int **range_id, size_t *no_range_id,
                  char **init_file, int **init_id, size_t *no_init_id,
                  char **basinhop_file, int *basinhop_id,
                  char **cubature_file, int *cubature_id, double **sim)
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
                // split multiple data files using :
                (*data_files) = strsplit(optarg, ":", no_files);
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
                    (*init_id) = strsplit_no(argv[optind], ":", no_init_id);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -i requires at least TWO arguments "
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
                else if (strcmp(optarg, "simplex") == 0)
                {
                	*method = 4;
                }
                else {
                    fprintf(stderr, "Invalid argument '%s': -o' option "
                            "requires one of these arguments: bfgs, "
                            "conj_pr, conj_fr, simplex \n\n",
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
                    (*range_id) = strsplit_no(argv[optind], ":", no_range_id);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -r requires at least TWO arguments "
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
                    (*group_id) = strsplit_no(argv[optind], ":", no_group_id);
                    optind++;
                }
                else
                {
                    fprintf(stderr, "Option -g requires at least TWO arguments "
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

    // verify and generate any errors regarding the arguments
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
    if (*model > 3 && (*init_id) == NULL)
    {
        fprintf(stderr,"Option -i is missing "
                "and it is required by -m D\n\n");
        req_flags[0] = 1;
    }

    // verify if init file is given when using -j
    if (*initial_values == TRUE && (*init_id) == NULL)
    {
        fprintf(stderr, "Option -i is missing "
                "and it is required by -j\n\n");
        req_flags[0] = 1;
    }

    // make sure initial_values is TRUE when init is given but -e is absent
    if ((*init_id) != NULL && *initial_estimation == FALSE)
    {
        *initial_values = TRUE;
    }

    // verify if number of data files matches with number of
    // init id, range id and group id
    if (*init_file != NULL && *no_init_id != 1 && *no_init_id != *no_files)
    {
        fprintf(stderr, "Number of IDs following -i has to "
                "either be 1 or the same as number of data files \n\n");
        req_flags[0] = 1;
    }
    if (*range_file != NULL && *no_range_id != 1 && *no_range_id != *no_files)
    {
        fprintf(stderr, "Number of IDs following -r has to "
                "either be 1 or the same as number of data files \n\n");
        req_flags[0] = 1;
    }
    if (*group_file != NULL && *no_group_id != 1 && *no_group_id != *no_files)
    {
        fprintf(stderr, "Number of IDs following -g has to "
                "either be 1 or the same as number of data files \n\n");
        req_flags[0] = 1;
    }

    if (req_flags[0] == 1)
    {
        print_usage(argv[0], stderr);
        return (EXIT_FAILURE);
    }

    // reallocate init id, range id and group id if needed
    int id = 0;
    if (*no_init_id == 1 && *no_files > 1)
    {
        id = (*init_id)[0];
        free(*init_id);
        *init_id = malloc(sizeof(int) * (*no_files));for (i = 0; i < (*no_files); i++)
        {
            (*init_id)[i] = id;
        }
    }
    if (*range_file != NULL && *no_range_id == 1 && *no_files > 1)
    {
        id = (*range_id)[0];
        free(*range_id);
        *range_id = malloc(sizeof(int) * (*no_files));
        for (i = 0; i < (*no_files); i++)
        {
            (*range_id)[i] = id;
        }
    }
    if (*group_file != NULL && *no_group_id == 1 && *no_files > 1)
    {
        id = (*group_id)[0];
        free(*group_id);
        *group_id = malloc(sizeof(int) * (*no_files));
        for (i = 0; i < (*no_files); i++)
        {
            (*group_id)[i] = id;
        }
    }

    if (*init_id == NULL)
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

int run_simulation(unsigned model, double *sim, char **data_files, size_t no_files,
             	   char *init_file, int *init_id, char *cubature_file, int cubature_id)
{
    int parsed = EXIT_SUCCESS;
    int status = EXIT_SUCCESS;
    char c;
    size_t i = 0;

    fprintf(stderr, "Warning: simulating data, ");
    for (i = 0; i < no_files; i++)
    {
        fprintf(stderr, "%s ", data_files[i]);
    }
    fprintf(stderr, "will be overwritten.\n");
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

    for (i = 0; i < no_files; i++)
    {
        // initialize the parameters
        ParamsModel pm;
        initialize_params_model(&pm);
        set_model(&pm, model);
        pm.n = sim[4];
        allocate_selection_params(&pm);
        allocate_grouping(&pm, NULL);

        parsed = parse(5, init_file, init_id[i], &pm, FALSE);
        status += handle_status(parsed, init_file, init_id[i]);

        if (status != EXIT_SUCCESS)
        {
            free_params_model(&pm);
            return (EXIT_FAILURE);
        }

        status = sim_data(&pm, sim, data_files[i]);
        free_params_model(&pm);

        if (status != EXIT_SUCCESS)
        {
            return (EXIT_FAILURE);
        }
    }

    return (EXIT_SUCCESS);
}

int run_optimization(unsigned model, unsigned kind, unsigned method,
             	 	 char **data_files, size_t no_files,
             	 	 char *group_file, int *group_id, char *optim_file, int optim_id,
             	 	 char *range_file, int *range_id, char *init_file, int *init_id,
             	 	 char *basinhop_file, int basinhop_id, char *cubature_file,
             	 	 int cubature_id, int verbose, int minutes,
					 int initial_estimation, int initial_values, int div)
{
    int status = EXIT_SUCCESS;
    int parsed = EXIT_SUCCESS;

    size_t i = 0;

    // parse optimization and bashin hopping parametrization
    ParamsOptim po;
    initialize_params_optim(&po);
    po.minutes = minutes;
    po.grad_descent = TRUE;

    ParamsBasinHop pb;
    initialize_params_basin_hop(&pb);

    parsed = parse(2, optim_file, optim_id, &po, initial_values);
    status += handle_status(parsed, optim_file, optim_id);

    parsed = parse(3, basinhop_file, basinhop_id, &pb, initial_values);
    status += handle_status(parsed, basinhop_file, basinhop_id);

    if (status != EXIT_SUCCESS)
    {
        return (EXIT_FAILURE);
    }

    // define optimization algorithm
    const void *type = NULL;
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
        case 4:
        {
        	type = gsl_multimin_fminimizer_nmsimplex2;
        	// I am not using a gradient descent
        	po.grad_descent = FALSE;
        	// update number of iterations - much larger for simplex
        	po.max_iter = 20000;
        	// update eps_abs - larger too
            po.eps_abs = 0.01;
			break;
        }
    }

    ParamsShare ps;
    allocate_params_share(&ps, no_files);
    ps.data_files = data_files;
    ps.verbose = verbose;
    ps.initial_estimation = initial_estimation;
    ps.initial_values = initial_values;
    if (init_id != NULL)
    {
        // if I don't use initial_estimation, make sure initial_values is TRUE
        if (ps.initial_estimation == FALSE)
        {
            ps.initial_values = TRUE;
        }
    }
    ps.kind = kind;

    printf("---- Performing inference using model ");
    switch (model)
    {
        case 1:
        {
            printf("A");
            break;
        }
        case 2:
        {
            printf("B");
            break;
        }
        case 3:
        {
            printf("C");
            break;
        }
        case 4:
        {
            printf("D");
            break;
        }
    }
    printf("\n\n");

    // parse data files and allocate params
    for (i = 0; i < no_files; i++)
    {
        printf("---- Performing inference on %s\n", data_files[i]);

        if (init_id != NULL)
        {
            printf("----      with -i %s %d\n", init_file, init_id[i]);
        }
        if (range_id != NULL)
        {
            printf("----      with -r %s %d\n", range_file, range_id[i]);
        }

        initialize_params(&ps.p[i]);
        set_model(ps.p[i].pm, model);
        // set various parameters that are already parsed
        ps.p[i].pm->div_flag = div;

        // the Params object is allocated in parse_data
        parsed = parse_data(data_files[i], &ps.p[i], &ps.p[i].pm->n);
        status += handle_status(parsed, data_files[i], 0);
        if (status != EXIT_SUCCESS)
        {
            break;
        }

        if (ps.p[i].counts_neut[0].sfs[ps.p[i].pm->n - 1] == -1)
        {
            printf("---- Data %s does not contain divergence counts. "
                   "Using option -w.\n",
                   data_files[i]);
            // make sure that div is not used
            ps.p[i].pm->div_flag = FALSE;
        }

        if (group_id != NULL)
        {
        	parsed = parse(6, group_file, group_id[i], ps.p[i].pm, initial_values);
        	status += handle_status(parsed, group_file, group_id[i]);
        }
        else
        {
        	// need to allocate default grouping
        	parsed = parse(6, group_file, -1, ps.p[i].pm, initial_values);
        	status += handle_status(parsed, group_file, -1);
        }

        if (range_file != NULL)
        {
            parsed = parse(0, range_file, range_id[i], ps.p[i].pm, initial_values);
            status += handle_status(parsed, range_file, range_id[i]);
        }

        if (status != EXIT_SUCCESS)
        {
            break;
        }

        initialize_selection_params(ps.p[i].pm, range_file);

        if (init_id != NULL)
        {
            parsed = parse(1, init_file, init_id[i], ps.p[i].pm, initial_values);
            status += handle_status(parsed, init_file, init_id[i]);
            // if failed because of wrong range and it's model D, add extra message
            if (parsed > WRONG_RANGE && ps.p[i].pm->model == 4)
            {
				fprintf(stderr,
						"Check that the given number of selection coefficients "
						"%d matches the one in %s with ID %d\n",
						ps.p[i].pm->no_sel / 2, init_file, init_id[i]);
            }
        }

        // check that I have at least two fragments for estimating a
        if (ps.p[i].no_neut == 1 && ps.p[i].no_sel == 1 && ps.p[i].pm->a != -1)
        {
            fprintf(stderr,
                    "---- Warning: mutation variability is not used for %s as "
                    "only one neutral and one selected fragment are available.\n",
                    data_files[i]);
            ps.p[i].pm->a = -1;
        }

        if (ps.p[i].pm->a == -1)
        {
            printf("---- No mutation variability. Using Poisson likelihood.\n");
            ps.p[i].pm->a_flag = FALSE;
            ps.p[i].pm->a_max = 1;
            ps.p[i].pm->a_min = -2;
        }

        // check range for the b parameter
        if (model < 4)
        {
        	if (ps.p[i].pm->sel_min[1] < 0.01)
        	{
        		fprintf(stderr,
						"---- Warning: the range for the b parameter for file %s\n"
						"              allows values below 0.01. When b is so "
        				"small, the DFE is very leptokurtic\n              "
        				"and it is very difficult to calculate the likelihood "
						"accurately.\n", data_files[i]);
        	}
        }

        if (status != EXIT_SUCCESS)
        {
            break;
        }

        printf("\n");
    }

    // start optimization
    if (status == EXIT_SUCCESS)
    {
        // first, update sharing
        set_sharing(&ps, init_file);

        // if r parameters are shared, make sure the grouping is the same
        if (ps.r_shared == TRUE)
        {
        	// for that, I need to make all possible comparisons
        	// and only check the groups that "overlap"
        	int min_r = 0;
        	int j = 0;
        	int k = 0;
        	for (i = 0; i < no_files - 1; i++)
        	{
        		for (j = i + 1; j < no_files; j++)
        		{
        			min_r = ps.p[i].pm->no_r > ps.p[j].pm->no_r
        					? ps.p[j].pm->no_r
							: ps.p[i].pm->no_r;
        			for (k = 0; k < min_r; k++)
        			{
        				if (ps.p[i].pm->inv_groups[k]
								   != ps.p[j].pm->inv_groups[k])
        				{
        					status = EXIT_FAILURE;
        					fprintf(stderr, "When the r parameters are shared and "
        							"option -g is used, the same groups should be "
        							"used for all files; files %s and %s have different groups",
									data_files[i], data_files[j]);
        					break;
        				}
        				if (status == EXIT_FAILURE)
        				{
        					break;
        				}
        			}
        		}
        	}
        }

        if (status == EXIT_SUCCESS)
        {
        	// make sure that if div is not used, the flag for lambda is FALSE
        	// and that if a == -1, the flag for a is FALSE
        	for (i = 0; i < no_files; i++)
        	{
        		if (ps.p[i].pm->div_flag == FALSE)
        		{
        			// do not estimate lambda
        			ps.p[i].pm->lambda_flag = FALSE;
        			ps.lambda_shared = FALSE;
        		}
        		if (ps.p[i].pm->a == -1)
        		{
        			// do not estimate a
        			ps.p[i].pm->a_flag = FALSE;
        			ps.a_shared = FALSE;
        		}
        	}

        	// count how many parameters I have to estimate
        	count_total_params(&ps);

        	allocate_params_basin_hop(&pb, ps);
        	run_basin_hopping(&pb, type, po, &ps);

        	// free memory
        	free_params_basin_hop(&pb);
        }
    }

    // free memory
    free_params_share(&ps);

    return (status);
}

void free_memory(double *sim, size_t no_files, char **data_files,
                 int *init_id, int *range_id, int *group_id)
{
    // size_t i = 0;

    if (sim)
    {
        free(sim);
    }

    if (init_id)
    {
        free(init_id);
    }

    if (range_id)
    {
        free(range_id);
    }

    if (group_id)
    {
        free(group_id);
    }
}

int main(int argc, char **argv)
{
    // this is for Eclipse working with stdout and stderr
    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    size_t i = 0;

    char **data_files = NULL;
    char *group_file = NULL;
    char *optim_file = NULL;
    char *range_file = NULL;
    char *init_file = NULL;
    char *basinhop_file = NULL;
    char *cubature_file = NULL;

    int *group_id = NULL;
    int optim_id = -1;
    int *range_id = NULL;
    int *init_id = NULL;
    int basinhop_id = -1;
    int cubature_id = -1;
    int initial_estimation = FALSE;
    int initial_values = FALSE;
    int div = TRUE; // use divergence data

    size_t no_files = 0;
    size_t no_group_id = 0;
    size_t no_range_id = 0;
    size_t no_init_id = 0;

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
                               &div, &data_files, &no_files,
                               &group_file, &group_id, &no_group_id,
							   &optim_file, &optim_id,
                               &range_file, &range_id, &no_range_id,
                               &init_file, &init_id, &no_init_id,
                               &basinhop_file, &basinhop_id, &cubature_file,
                               &cubature_id, &sim);

    if (status != EXIT_SUCCESS)
    {
        free_memory(sim, no_files, data_files, init_id, range_id, group_id);

        return(status);
    }

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
        status = run_simulation(model, sim, data_files, no_files, init_file,
        		                init_id, cubature_file, cubature_id);
    }
    else
    {
        status = run_optimization(model, kind, method, data_files, no_files,
        		                  group_file, group_id, optim_file, optim_id,
								  range_file, range_id, init_file, init_id,
								  basinhop_file, basinhop_id, cubature_file,
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

    free_memory(sim, no_files, data_files, init_id, range_id, group_id);

    return (status);
}
