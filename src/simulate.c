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

#include "simulate.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>
#include <string.h>

#include "likelihood.h"

void initialize_random_generator(gsl_rng **rand)
{
    // random number generator
    const gsl_rng_type * T;
    gsl_rng_env_setup();

    // set the seed
    struct timeval tv;
    gettimeofday(&tv, 0);
    gsl_rng_default_seed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default;
    *rand = gsl_rng_alloc(T);
}

void sim_frag(FILE *f, ParamsModel pm, double len, double *expec, gsl_rng *rand)
{
    double p = 0;
    unsigned i = 0;
    int aux = 0;

    for (i = 0; i < pm.n - 1; i++)
    {
        p = pm.a / (pm.a + expec[i] * len);
        aux = gsl_ran_negative_binomial(rand, p, pm.a);
        fprintf(f, "%d\t", aux);
    }
    fprintf(f, "%g\t\t", len);

    // divergence counts need to be written at the end
    i = pm.n - 1;
    p = pm.a / (pm.a + expec[i] * len);
    aux = gsl_ran_negative_binomial(rand, p, pm.a);
    fprintf(f, "%d\t", aux);
    fprintf(f, "%g\n", len);
}

int sim_data(ParamsModel *pm, double *len, char *filename)
{
    FILE *f = fopen(filename, "w");
    if (f == NULL)
    {
        fprintf(stderr, "Simulation failed due to fopen on %s: %s\n",
                filename, strerror(errno));
        return (EXIT_FAILURE);
    }

    fprintf_params_model(*pm, f, "#");

    unsigned i = 0;

    // initialize arrays for storing expectations
    double *sel = malloc(sizeof(double) * pm->n);
    double *neut = malloc(sizeof(double) * pm->n);
    set_to_zero(&sel, pm->n);
    set_to_zero(&neut, pm->n);

    set_sel_expec(pm, &sel, FALSE);
    set_neut_expec(*pm, &neut);

    // add ancestral error
    set_anc_expec(pm->eps_an, pm->n, &sel);
    set_anc_expec(pm->eps_an, pm->n, &neut);

    // random number generator
    gsl_rng *rand = NULL;
    initialize_random_generator(&rand);

    fprintf(f, "\n%g %g %d\n\n", len[0], len[2], pm->n);

    // simulate neutral fragments
    for (i = 0; i < len[0]; i++)
    {
        sim_frag(f, *pm, len[1], neut, rand);
    }
    fprintf(f, "\n");

    // simulate selected fragments
    for (i = 0; i < len[2]; i++)
    {
        sim_frag(f, *pm, len[3], sel, rand);
    }

    // free memory
    free(sel);
    free(neut);
    gsl_rng_free(rand);

    fclose(f);

    return (EXIT_SUCCESS);
}
