// S2LET package
// Copyright (C) 2012-2014
// Boris Leistedt & Jason McEwen

#ifndef S2LET_TYPES
#define S2LET_TYPES

#include "ssht.h"

typedef enum {
    /*!
     * McEwen and Wiaux sampling:
     * 2*L-1 samples in alpha, in [0, 2pi).
     * L samples in beta, in (0, pi].
     * 2*N-1 samples in gamma, in [0, 2pi).
     */
    S2LET_SAMPLING_MW,
    /*!
     * McEwen and Wiaux symmetric sampling:
     * 2*L samples in alpha, in [0, 2pi).
     * L+1 samples in beta, in [0, pi].
     * 2*N-1 samples in gamma, in [0, 2pi).
     */
    S2LET_SAMPLING_MW_SS,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    S2LET_SAMPLING_SIZE
} s2let_sampling_t;

typedef enum {
    /*! Use the standard normalisation of wavelets such that they fulfil the
     *  the identity relation */
    S2LET_WAV_NORM_DEFAULT,
    /*! Normalise the wavelets as if they were spin-lowered. */
    S2LET_WAV_NORM_SPIN_LOWERED,
    /*!
     * "guard" value that equals the number of usable enum values.
     * useful in loops, for instance.
     */
    S2LET_WAV_NORM_SIZE
} s2let_wav_norm_t;


/*!
 * A struct with all parameters that are common to several
 * functions of the API. In general only one struct needs to
 * be created and a const pointer to it is passed around.
 * \attention
 *   Make sure to use an initializer upon
 *   declaration, even if it is left empty. This
 *   ensures that all fields are initialized to
 *   zero (even if in non-static storage). This way,
 *   your code will remain backwards compatible if
 *   more fields are added to this struct in the
 *   future:
 *   \code{.c}
 *   s2let_parameters_t parameters = {};
 *   \endcode
 */
typedef struct {
    /*!
     * Detail level for diagnostic console output in range [0,5].
     * \var int verbosity
     */
    int verbosity;

    /*!
     * A non-zero value indicates that the signal f is real. Not
     * all functions respect this value - instead there may be
     * separate complex and real functions. See the documentation
     * of each function for details.
     * \var int reality
     */
    int reality;

    /*!
     * A non-zero value indicates that the signal is stored in a
     * full-resolution format, where each wavelet scale is
     * upsampled to use the same amount of pixels.
     * A zero-value indicates that the signal uses only as many
     * pixels as necessary for each wavelet scale's upper harmonic
     * bandlimit. This can lead to significant storage and time
     * savings and is the default behaviour.
     * \var int upsample
     */
    int upsample;

    /*!
     * Wavelet parameter which determines the scale factor between
     * consecutive wavelet scales.
     * \var int B
     */
    int B;

    /*!
     * Upper harmonic band-limit. Only flmn with l < L will be stored
     * and considered.
     * \var int L
     */
    int L;

    /*!
     * First wavelet scale to be used.
     * \var int J_min
     */
    int J_min;

    /*!
     * Upper azimuthal band-limit. Only flmn with n < N will
     * be stored.
     * \var int N
     */
    int N;

    /*!
     * Spin number of the signal f.
     * \var int spin
     */
    int spin;

    /*!
     * Indicates how to normalise the wavelets and scaling function.
     * \var s2let_wav_norm_t normalization
     */
    s2let_wav_norm_t normalization;
    


    /*!
     * If normalization has value S2LET_WAV_NORM_SPIN_LOWERED
     * this parameter indicates which spin number the wavelets
     * were lowered from. Otherwise, it is ignored.
     */
    int original_spin;

    /*!
     * Sampling scheme to use for samples of the signal f as well as
     * the wavelets.
     * \var s2let_sampling_t sampling_scheme
     */
    s2let_sampling_t sampling_scheme;

    /*!
     * Recursion method to use for computing Wigner functions.
     * \var ssht_dl_method_t dl_method
     */
    ssht_dl_method_t dl_method;

} s2let_parameters_t;

#endif
