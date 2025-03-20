/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.s2let;

import org.junit.Test;
import static org.junit.Assert.*;

/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */
public class AxisymmetricWaveletTransformTest {
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_harm() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_harm(32, 2, 1, false, false, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(32, 3, 1, true, false, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(64, 2, 2, false, true, verbosity);
        test_AxisymmetricWaveletTransformTest_harm(32, 3, 0, true, true, verbosity);
    }
    /**
     * Test the exactness of the axisymmetric wavelet transform in harmonic space.
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param multires the multiresolution flag, i.e. <code>true</code> if all wavelet scales need to be constructed at the resolution corresponding to their bandlimit
     * @param reality the reality flag, i.e. weither or not the harmonic coefficient correspond to real signals on the sphere
     * @param verbosity 
     */
    public void test_AxisymmetricWaveletTransformTest_harm(int bandlimit, int waveletParameter, int firstScale, boolean multires, boolean reality, boolean verbosity) {
        SphericalHarmonicTransform flm;
        if (!reality) {
            flm = SphericalHarmonicTransform.random(bandlimit);
        } else {
            flm = SphericalHarmonicTransform.randomReal(bandlimit);
        }
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromHarmonics(flm, waveletParameter, firstScale, multires); 
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromWavelets(wavelets);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of Harm-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", multires="+multires+", reality="+reality+" : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_MWMap() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_MWMap(32, 2, 1, false, false, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(64, 3, 1, true, false, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(128, 2, 3, false, true, verbosity);
        test_AxisymmetricWaveletTransformTest_MWMap(32, 2, 0, true, true, verbosity);
    }
    /**
     * Test the exactness of the axisymmetric wavelet transform in real space: MW sampling.
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param multires the multiresolution flag, i.e. <code>true</code> if all wavelet scales need to be constructed at the resolution corresponding to their bandlimit
     * @param reality the reality flag, i.e. weither or not the harmonic coefficient correspond to real signals on the sphere
     * @param verbosity 
     */
    public void test_AxisymmetricWaveletTransformTest_MWMap(int bandlimit, int waveletParameter, int firstScale, boolean multires, boolean reality, boolean verbosity) {
        SphericalHarmonicTransform flm;
        if (!reality) {
            flm = SphericalHarmonicTransform.random(bandlimit);
        } else {
            flm = SphericalHarmonicTransform.randomReal(bandlimit);
        }        
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.MW);
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromMWMap(map, waveletParameter, firstScale, multires); 
        PixelizedMap map_rec = PixelizedMap.fromWavelets(wavelets, SamplingScheme.MW);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map_rec);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of MW-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", multires="+multires+", reality="+reality+" : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_AxisymmetricWaveletTransformTest_HealpixMap() {
        boolean verbosity = true;
        test_AxisymmetricWaveletTransformTest_HealpixMap(32, 2, 1, 64, verbosity);
        test_AxisymmetricWaveletTransformTest_HealpixMap(128, 2, 3, 128, verbosity);
    }
    
    /**
     * Test that the axisymmetric wavelet transform works for Healpix maps.
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param nside the resolution of the Healpix maps
     * @param reality the reality flag, i.e. weither or not the harmonic coefficient correspond to real signals on the sphere
     * @param verbosity 
     */
    public void test_AxisymmetricWaveletTransformTest_HealpixMap(int bandlimit, int waveletParameter, int firstScale, int nside, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.randomReal(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.HEALPIX, nside);
        AxisymmetricWaveletTransform wavelets = AxisymmetricWaveletTransform.fromHealpixMap(map, waveletParameter, firstScale, bandlimit); 
        PixelizedMap map_rec = PixelizedMap.fromWavelets(wavelets, SamplingScheme.HEALPIX, nside);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map_rec, bandlimit);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if (verbosity) {
            System.out.println("MaxAbsErr of HPX-AxiWav with L=" + bandlimit + ", B="+waveletParameter+", J_min="+firstScale+", nside="+nside+" : " + error);
        }
        assertEquals(error, 0.0, 1.0);
    }
    
}
