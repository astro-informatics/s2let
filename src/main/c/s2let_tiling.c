// S2LET package
// Copyright (C) 2012
// Boris Leistedt & Jason McEwen

#include "s2let.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//
//
typedef enum { S2DW, NEEDLET, SPLINE } s2let_kernel_type ;
s2let_kernel_type s2let_kernel = S2DW;
//
//

/*!
 * Switch to different wavelet type.
 *
 * \param[in]  typenum Integer: 1 for scale-discretised, 2 for needlets and 3 for spline wavelets.
 * \retval none
 */
void s2let_switch_wavtype(int typenum)
{
    //printf("Input wavelet type : %i\n",typenum);
    if (typenum == 1){
    //printf("Kernel switch 1: using scale-discretised wavelets.\n");
        s2let_kernel = S2DW;
    } else if (typenum == 2){
    //printf("Kernel switch 2: using needlets.\n");
        s2let_kernel = NEEDLET;
    } else if (typenum == 3){
    //printf("Kernel switch 3: using cubic splines wavelets.\n");
        s2let_kernel = SPLINE;
    } else {
        printf("Kernel number should be 1, 2 or 3. Default is 1.\n");
        s2let_kernel = S2DW;
    }
}

/*!
 * Computes band-limit of a specific wavelet scale.
 *
 * \param[in]  j Wavelet scale.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 * \retval band-limit
 */
int s2let_bandlimit(int j, const s2let_parameters_t *parameters)
{
    double B = parameters->B;
    int L = parameters->L;
    // int J_min = parameters->J_min;

    int Jmax;
    switch (s2let_kernel)
    {
    case S2DW:
    case NEEDLET:
        return ceil(pow(B, j+1));
    case SPLINE:
        Jmax = s2let_j_max(parameters);
        if (j == Jmax) return L;
        //if (j < J_min) return ceil(L / (double) pow(B, Jmax-J_min-1));
        return ceil(L / pow(B, Jmax-j-2));
    default:
        // This should never happen
        return -1;
    }
}

/*!
 * Computes the minimum harmonic index supported by the given
 * wavelet scale.
 *
 * \param[in]  j Wavelet scale.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink
 * \retval el_min
 */
int s2let_L0(int j, const s2let_parameters_t *parameters)
{
    double B = parameters->B;

    switch (s2let_kernel)
    {
    case S2DW:
    case NEEDLET:
        return ceil(pow(B, j-1));
    case SPLINE:
        return 0;
    default:
        // This should never happen
        return -1;
    }
}

/*!
 * Computes needlet maximum level required to ensure exact reconstruction.
 *
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink
 * \retval j_max
 */
int s2let_j_max(const s2let_parameters_t *parameters)
{
  double B = parameters->B;
  int L = parameters->L;

  return ceil(log(L) / log(B));
}

/*!
 * Allocates axisymmetric tiling kernels in harmonic space.
 *
 * \param[out]  kappa Kernel functions for the wavelets.
 * \param[out]  kappa0 Kernel for the scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink
 * \retval none
 */
void s2let_tiling_axisym_allocate(double **kappa, double **kappa0, const s2let_parameters_t *parameters)
{
    int L = parameters->L;

    int J = s2let_j_max(parameters);
    *kappa = calloc((J+1) * L, sizeof **kappa);
    *kappa0 = calloc(L, sizeof **kappa0);
}

void s2let_tiling_phi2_s2dw(double *phi2, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    double B = parameters->B;

    int j, l;
    int J = s2let_j_max(parameters);
    int n = 300;

    double kappanorm = s2let_math_kappa0_quadtrap_s2dw(1.0 / B, 1.0, n, B);
    for (j = 0; j <= J+1; j++){
        for (l = 0; l < L; l++){
            if (l < pow(B,j-1)) {
                phi2[l+j*L] = 1;
            } else if (l > pow(B,j)) {
                phi2[l+j*L] = 0;
            } else {
                phi2[l+j*L] = s2let_math_kappa0_quadtrap_s2dw((double)l / pow(B, j), 1.0, n, B) / kappanorm;
            }
        }
    }
}

void s2let_tiling_phi2_needlet(double *phi2, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    double B = parameters->B;

    int j, l;
    int J = s2let_j_max(parameters);
    int n = 300;
    double u;

    double kappanorm = s2let_math_kappa0_quadtrap_needlet(-1.0, 1.0, n);
    for (j = 0; j <= J+1; j++){
        for (l = 0; l < L; l++){
            if (l < pow(B,j-1)) {
                phi2[l+j*L] = 1;
            } else if (l > pow(B,j)) {
                phi2[l+j*L] = 0;
            } else {
                u = 1.0 - 2.0 * B / (B - 1.0) * ( l * pow(B, -j) - 1.0 / B );
                phi2[l+j*L] = s2let_math_kappa0_quadtrap_needlet(-1.0, u, n) / kappanorm;
            }
        }
    }
}

void s2let_tiling_phi2_spline(double *phi2, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    double B = parameters->B;

    int j = 0, l;
    int J = s2let_j_max(parameters);
    phi2[(J+1-j)*L] = 1.0;
    for (l = 1; l < L; l++){
        phi2[l+(J+1-j)*L] = 1.0;
    }
    for (j = 1; j <= J+1; j++){
        double bl = (double) L / (double) pow(B, j-2);
        phi2[(J+1-j)*L] = 1.0;
        for (l = 1; l < L; l++){
            if (l > bl)
                phi2[l+(J+1-j)*L] = 0.0;
            else
                phi2[l+(J+1-j)*L] = s2let_math_spline_scalingfct((double) l, bl);
        }
    }
}

/*!
 * Generates axisymmetric tiling in harmonic space.
 *
 * \param[out]  kappa Kernel functions for the wavelets.
 * \param[out]  kappa0 Kernel for the scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 * \retval none
 */
void s2let_tiling_axisym(double *kappa, double *kappa0, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    int J_min = parameters->J_min;

    int j, l;
    int J = s2let_j_max(parameters);

    double previoustemp = 0.0, temp;
    double *phi2 = (double*)calloc((J+2) * L, sizeof(double));

    if(s2let_kernel == SPLINE)
        s2let_tiling_phi2_spline(phi2, parameters); // SPLINE tiling
    if(s2let_kernel == S2DW)
        s2let_tiling_phi2_s2dw(phi2, parameters); // S2DW tiling
    if(s2let_kernel == NEEDLET)
        s2let_tiling_phi2_needlet(phi2, parameters); // Needlet tiling

    for (l = 0; l < L; l++)
        kappa0[l] = sqrt(phi2[l+J_min*L]);

    for (j = J_min; j <= J; j++){
        for (l = 0; l < L; l++){
            temp = sqrt(phi2[l+(j+1)*L] - phi2[l+j*L]);
            if( isnan(temp) || isinf(temp) )
                kappa[l+j*L] = previoustemp;
            else
                kappa[l+j*L] = temp;
            previoustemp = temp;
        }
        for (l = 0; l < L; l++)
            if( !isfinite(kappa[l+j*L]) )
                kappa[l+j*L] = kappa[l+j*L-1];
    }
    free(phi2);
}

/*!
 * Allocates space for directionality components in harmonic
 * space.
 *
 * \param[out]  s_elm Pointer to allocated space for harmonic
 *                    coefficients of directionality components.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_tiling_direction_allocate(complex double **s_elm, const s2let_parameters_t *parameters)
{
  int L = parameters->L;
  // TODO: This could be reduced by not storing s_elm with |m| >= N
  *s_elm = calloc(L*L, sizeof **s_elm);
}

/*!
 * Generates the harmonic coefficients for the directionality
 * component of the tiling functions.
 * This implementation is based on equation (11) in the wavelet
 * computation paper.
 *
 * \param[out]  s_elm Harmonic coefficienets of directionality
 *                    components.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::N N\endlink
 *                        \link s2let_parameters_t::spin spin\endlink
 *
 */
void s2let_tiling_direction(complex double *s_elm, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    int N = parameters->N;

    // TODO: Add spin parameter to avoid computation of el < |s|
    complex double nu;
    int el, m, ind;

    if (N % 2)
        nu = 1;
    else
        nu = I;

    // Skip the s_00 component, as it is zero.
    ind = 1;

    for (el = 1; el < L; ++el)
    {
        int gamma;
        // This if else replaces the -1^(N+l)
        if ((N+el) % 2)
            gamma = MIN(N-1, el);
        else
            gamma = MIN(N-1, el-1);

        for (m = -el; m <= el; ++m)
        {
            // This if/else takes care of the azimuthal
            // band-limit and replaces the beta factor.
            if (ABS(m) < N && (N+m) % 2)
                s_elm[ind] = nu*sqrt(binomial_coefficient(gamma, (gamma - m)/2UL, 1)/pow(2, gamma));
            else
                s_elm[ind] = 0.0;

            ++ind;
        }
    }
}

/*!
 * Allocates space for directional wavelets in harmonic space.
 *
 * \param[out]  psi Pointer to allocated space for harmonic
 *                  coefficients of directional wavelets.
 * \param[out]  phi Pointer to allocated space for harmonic
 *                  coefficients of scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::N N\endlink
 * \retval none
 */
void s2let_tiling_wavelet_allocate(complex double **psi, double **phi, const s2let_parameters_t *parameters)
{
    int L = parameters->L;

    // TODO: This could be reduced by not storing psi_j_elm with |m| >= N
    int J = s2let_j_max(parameters);
    *psi = calloc((J+1) * L*L, sizeof **psi);
    *phi = calloc(L, sizeof **phi);
}

/*!
 * Computes the normalization factor for spin-lowered wavelets,
 * which is sqrt((l+s)!/(l-s)!).
 *
 * \param[in]  el    Harmonic index el.
 * \param[in]  spin  Spin number the wavelet was lowered from.
 */
static double s2let_spin_normalization(int el, int spin)
{
    double factor = 1;
    int s;

    for (s = -ABS(spin)+1; s <= ABS(spin); ++s)
    {
        factor *= el+s;
    }

    if (spin > 0)
        return sqrt(factor);
    else
        return sqrt(1.0/factor);
}

/*!
 * Generates the harmonic coefficients for the directional tiling wavelets.
 * This implementation is based on equation (7) in the wavelet
 * computation paper.
 *
 * \param[out]  psi Harmonic coefficienets of directional wavelets.
 * \param[out]  phi Harmonic coefficienets of scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 *                        \link s2let_parameters_t::spin spin\endlink
 *                        \link s2let_parameters_t::original_spin original_spin\endlink
 *
 */
void s2let_tiling_wavelet(complex double *psi, double *phi, const s2let_parameters_t *parameters) {
    int L = parameters->L;
    int J_min = parameters->J_min;
    int spin = parameters->spin;
    int original_spin = parameters->original_spin;

    // TODO: Add spin parameter to avoid computation of el < |s|
    // TODO: Correctly compute spin scaling functions
    double *kappa;
    double *kappa0;
    complex double *s_elm;
    int j, el, m, el_min;
    int J = s2let_j_max(parameters);

    // TODO: Allocate kappa0 directly inside phi. For this, we should probably
    //       separate the allocation functions to do only one allocation per
    //       function.
    s2let_tiling_axisym_allocate(&kappa, &kappa0, parameters);
    s2let_tiling_axisym(kappa, kappa0, parameters);
    s2let_tiling_direction_allocate(&s_elm, parameters);
    s2let_tiling_direction(s_elm, parameters);

    el_min = MAX(ABS(spin), ABS(original_spin));

    for (el = el_min; el < L; ++el)
    {
        phi[el] = sqrt((2*el+1)/(4.0*PI)) * kappa0[el];
        if (original_spin != 0)
            phi[el] *= s2let_spin_normalization(el, original_spin) * pow(-1, original_spin);
    }

    for (j = J_min; j <= J; ++j)
    {
        int ind = el_min*el_min;
        for (el = el_min; el < L; ++el)
        {
            for (m = -el; m <= el; ++m)
            {
                psi[j*L*L + ind] = sqrt((2*el+1)/(8.0*PI*PI)) * kappa[j*L + el] * s_elm[ind];
                if (original_spin != 0)
                    psi[j*L*L + ind] *= s2let_spin_normalization(el, original_spin) * pow(-1, original_spin);
                ++ind;
            }
        }
    }

    free(kappa);
    free(kappa0);
    free(s_elm);
}

/*!
 * Checks exactness of the harmonic tiling kernels by checking
 * the admissibility condition.
 *
 * \param[in]  kappa Kernel functions for the wavelets.
 * \param[in]  kappa0 Kernel for the scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *
 * \retval Achieved accuracy (should be lower than e-14).
 */
double s2let_tiling_axisym_check_identity(double *kappa, double *kappa0, const s2let_parameters_t *parameters)
{
    int L = parameters->L;

    int l, j;
    int J = s2let_j_max(parameters);
    //int l_min = s2let_el_min(B, J_min);
    double error = 0;

    double *ident;
    ident = calloc(L, sizeof *ident);

    for (l = 0; l < L; l++)
        ident[l] = pow(kappa0[l], 2.0);

    for (l = 0; l < L; l++) {
        for (j = 0; j <= J; j++) {
            ident[l] += pow(kappa[l+j*L], 2.0);
        }

        error = MAX(error, fabs(ident[l] - 1.0));
    }

    free(ident);
    return error;
}

/*!
 * Checks exactness of the directionality components by checking
 * the admissibility condition.
 *
 * \param[in]  s_elm Harmonic coefficients of directionality
 *                   components.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::N N\endlink
 * \retval Achieved accuracy (should be lower than e-14).
 */
double s2let_tiling_direction_check_identity(complex double *s_elm, const s2let_parameters_t *parameters)
{
    int L = parameters->L;

    int el, m, ind;
    double error = 0.0; // maximum error for all el

    // Skip the s_00 component, as it is zero.
    ind = 1;

    for (el = 1; el < L; ++el)
    {
        double sum = 0.0; // sum for each el
        for (m = -el; m <= el; ++m)
        {
            sum += s_elm[ind] * conj(s_elm[ind]);
            ++ind;
        }

        error = MAX(error, fabs(sum - 1.0));
    }

    return error;
}

/*!
 * Checks exactness of the directional wavelets by checking
 * the admissibility condition.
 *
 * \param[in]  psi Harmonic coefficients of directional wavelets.
 * \param[in]  phi Harmonic coefficients of scaling function.
 * \param[in]  parameters A parameters object with (at least) the following fields:
 *                        \link s2let_parameters_t::B B\endlink,
 *                        \link s2let_parameters_t::L L\endlink,
 *                        \link s2let_parameters_t::J_min J_min\endlink
 *                        \link s2let_parameters_t::N N\endlink
 *                        \link s2let_parameters_t::spin spin\endlink
 * \retval Achieved accuracy (should be lower than e-14).
 */
double s2let_tiling_wavelet_check_identity(complex double *psi, double *phi, const s2let_parameters_t *parameters)
{
    int L = parameters->L;
    int spin = parameters->spin;

    int j, el, m, ind;
    int J = s2let_j_max(parameters);
    double error = 0.0; // maximum error for all el

    double *ident;
    ident = calloc(L, sizeof *ident);

    for (el = ABS(spin); el < L; ++el)
    {
        ident[el] += 4.0*PI/(2*el+1) * phi[el] * phi[el];
    }

    for (j = 0; j <= J; ++j)
    {
        ind = spin*spin;
        for (el = ABS(spin); el < L; ++el)
        {
            for (m = -el; m <= el; ++m)
            {
                ident[el] += 8.0*PI*PI/(2*el+1) *
                             psi[j*L*L + ind] * conj(psi[j*L*L + ind]);
                ++ind;
            }
        }
    }

    for (el = ABS(spin); el < L; ++el)
    {
        error = MAX(error, fabs(ident[el] - 1.0));
    }

    return error;
}
