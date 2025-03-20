
package org.s2let;

import org.junit.*;
import static org.junit.Assert.*;


/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */
public class MWMapTest {
    
    @Test
    public void test_MW_SphericalHarmonicTransform() {
        boolean verbosity = true;
        test_MW_SphericalHarmonicTransform(16, verbosity);
        test_MW_SphericalHarmonicTransform(64, verbosity);
        test_MW_SphericalHarmonicTransform(128, verbosity);
    }
    /**
     * Test the exactness of the spherical harmonic transform built on the MW sampling (complex signals).
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param verbosity 
     */
    public void test_MW_SphericalHarmonicTransform(int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.random(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.MW);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of MW-SHA with L=" + bandlimit + " : " + error); 
        }
        assertEquals(error, 0.0, 1e-8);
    }
    
    @Test
    public void test_MWReal_SphericalHarmonicTransform() {
        boolean verbosity = true;
        test_MWReal_SphericalHarmonicTransform(16, verbosity);
        test_MWReal_SphericalHarmonicTransform(64, verbosity);
        test_MWReal_SphericalHarmonicTransform(128, verbosity);
    }
    /**
     * Test the exactness of the spherical harmonic transform built on the MW sampling (real signals).
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param verbosity 
     */
    public void test_MWReal_SphericalHarmonicTransform(int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.randomReal(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.MW);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of RealMW-SHA with L=" + bandlimit + " : " + error);
        }
        assertEquals(error, 0.0, 1e-8);
    }
}
