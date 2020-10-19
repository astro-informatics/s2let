// S2LET package
// Copyright (C) 2014
// Boris Leistedt & Jason McEwen

#ifndef S2LET_SO3
#define S2LET_SO3

#include <so3/so3.h>
#include <time.h>

// Define a few macros for fixed SO3 configuration used
// throughout S2LET.
#define S2LET_SO3_N_ORDER SO3_N_ORDER_NEGATIVE_FIRST
#define S2LET_SO3_STORAGE SO3_STORAGE_COMPACT

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * A static helper function to prepopulate an so3_parameters_t
 * struct with data from an s2let_parameters_t struct.
 */
static inline void fill_so3_parameters(so3_parameters_t *so3_parameters, const s2let_parameters_t *parameters)
{
    so3_parameters->verbosity = parameters->verbosity;
    so3_parameters->L = parameters->L;
    so3_parameters->N = parameters->N;
    so3_parameters->sampling_scheme = (so3_sampling_t) parameters->sampling_scheme;
    so3_parameters->n_order = S2LET_SO3_N_ORDER;
    so3_parameters->storage = S2LET_SO3_STORAGE;
    so3_parameters->dl_method = parameters->dl_method;
    so3_parameters->reality = parameters->reality;
    so3_parameters->steerable = 0;  // When direct routines support steerable so3 change this to 1.
    so3_parameters->L0 = 0;

    if (parameters->N % 2)
        so3_parameters->n_mode = SO3_N_MODE_EVEN;
    else
        so3_parameters->n_mode = SO3_N_MODE_ODD;
    //so3_parameters->n_mode = SO3_N_MODE_ALL; // just to test
}

static inline void print_so3_parameters(so3_parameters_t *so3_parameters)
{
    FILE *f;
    f = fopen("test.txt", "a");

    fprintf(f, "hello, %f\n", (double)clock());
    fclose(f);

    fprintf (stdout,"so3_parameters->verbosity        %d\n", so3_parameters->verbosity);
    fprintf (stdout,"so3_parameters->reality          %d\n", so3_parameters->reality);
    fprintf (stdout,"so3_parameters->L0               %d\n", so3_parameters->L0);
    fprintf (stdout,"so3_parameters->L                %d\n", so3_parameters->L);
    fprintf (stdout,"so3_parameters->N                %d\n", so3_parameters->N);
    fprintf (stdout,"so3_parameters->sampling_scheme  %d\n", (int)so3_parameters->sampling_scheme);
    fprintf (stdout,"so3_parameters->n_order          %d\n", (int)so3_parameters->n_order);
    fprintf (stdout,"so3_parameters->storage          %d\n", (int)so3_parameters->storage);
    fprintf (stdout,"so3_parameters->n_mode           %d\n", (int)so3_parameters->n_mode);
    fprintf (stdout,"so3_parameters->dl_method        %d\n", (int)so3_parameters->dl_method);
    fprintf (stdout,"so3_parameters->steerable        %d\n", so3_parameters->steerable);
    fflush (stdout);
}

#ifdef __cplusplus
}
#endif
#endif
