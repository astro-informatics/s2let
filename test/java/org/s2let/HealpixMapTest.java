
package org.s2let;

import org.bridj.Pointer;
import org.junit.*;
import static org.junit.Assert.*;

/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */
public class HealpixMapTest {
    
    @Test
    public void test_Healpix_SphericalHarmonicTransform() {
        boolean verbosity = true;
        test_Healpix_SphericalHarmonicTransform(64, 16, verbosity);
        test_Healpix_SphericalHarmonicTransform(256, 128, verbosity);
    }
    /**
     * Test that the spherical harmonic transform works well on the Healpix sampling (real signals).
     * @param nside the resolution of the Healpix maps
     * @param bandlimit the bandlimit of the decomposition to generate and work with
     * @param verbosity 
     */
    public void test_Healpix_SphericalHarmonicTransform(int nside, int bandlimit, boolean verbosity) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.randomReal(bandlimit);
        PixelizedMap map = PixelizedMap.fromHarmonics(flm, SamplingScheme.HEALPIX, nside);
        SphericalHarmonicTransform flm_rec = SphericalHarmonicTransform.fromMap(map, bandlimit);
        double error = flm_rec.maxAbsoluteDifferenceWith(flm);
        if(verbosity){
            System.out.println("MaxAbsErr of HPX-SHA with L=" + bandlimit + " and nside=" +nside + " : " + error);
        }
        assertEquals(error, 0.0, 1e-1);
    }
    
}
