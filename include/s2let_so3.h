// S2LET package
// Copyright (C) 2014
// Boris Leistedt & Jason McEwen

#ifndef S2LET_SO3
#define S2LET_SO3

#include <so3.h>

// Define a few macros for fixed SO3 configuration used
// throughout S2LET.
#define S2LET_SO3_N_ORDER SO3_N_ORDER_NEGATIVE_FIRST
#define S2LET_SO3_STORAGE SO3_STORAGE_COMPACT
// SO3_STORAGE_PADDED   SO3_STORAGE_COMPACT

/*!
 * A static helper function to prepopulate an so3_parameters_t
 * struct with data from an s2let_parameters_t struct.
 */
static inline void fill_so3_parameters(so3_parameters_t *so3_parameters, const s2let_parameters_t *parameters)
{
    so3_parameters->verbosity = parameters->verbosity;
    so3_parameters->L = parameters->L;
    so3_parameters->N = parameters->N;
    so3_parameters->sampling_scheme = parameters->sampling_scheme;
    so3_parameters->n_order = S2LET_SO3_N_ORDER;
    so3_parameters->storage = S2LET_SO3_STORAGE;
    so3_parameters->dl_method = parameters->dl_method;
    so3_parameters->reality = parameters->reality;
    so3_parameters->steerable = 1;

    if (parameters->N % 2)
        so3_parameters->n_mode = SO3_N_MODE_EVEN;
    else
        so3_parameters->n_mode = SO3_N_MODE_ODD;
}

#endif
