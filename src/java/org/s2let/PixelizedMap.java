
package org.s2let;

import org.bridj.*;
import static org.s2let.bindings.S2letLibrary.*;

/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */

/**
 * Abstract class collecting the common properties of the Healpix and MW samplings
 * and providing public interfaces / constructors, so that the user can avoid to 
 * manipulate different classes explicitly.
 */
public abstract class PixelizedMap {
    /**
     * Code name for the sampling scheme.
     */
    protected final SamplingScheme scheme;
    /**
     * Resolution of the map: nside for Healpix, L for MW (which is also the bandlimit).
     */
    protected final int resolution;
    /**
     * Reality flag: if the map is real, then it contains doubles instead of complex doubles.
     */
    protected final boolean reality;
    /**
     * Getter for the scheme code name.
     * @return the scheme flag
     */
    public SamplingScheme getScheme() {
        return scheme;
    }
    /**
     * Getter for the reality flag.
     * @return <code>true</code> if the signal is real, <code>false</code> if it is complex
     */
    public boolean getReality() {
        return this.reality;
    };
    /**
     * Getter for the resolution field.
     * @return the resolution
     */
    public int getResolution() {
        return this.resolution;
    };
    /**
     * Private constructor to initialize the simple fields.
     * @param sc the flag of the sampling scheme
     * @param reso the resolution of the map
     * @param rea the reality flag
     */
    protected PixelizedMap(SamplingScheme sc, int reso, boolean rea) {
        this.scheme = sc;
        this.reality = rea;
        this.resolution = reso;
    }
    /**
     * Factory method to create a Pixelized map from a spherical harmonic decomposition, a sampling scheme flag and the resolution.
     * @param flm the spherical harmonic decomposition
     * @param sc the flag of the sampling scheme
     * @param nside the Healpix resolution, only required if the map is Healpix
     * @return a new instance for either a MW or a Healpix map
     */
    static public PixelizedMap fromHarmonics(SphericalHarmonicTransform flm, SamplingScheme sc, int nside) {
        if (sc == SamplingScheme.MW){        
            if (flm.getReality()) {
                return new MWRealMap(flm);
            } else {
                return new MWCmplxMap(flm);
            }
        } else {
            return new HealpixRealMap(flm, nside);
        }
            
    };
    /**
     * Simplified factory method to create a Pixelized map from a spherical harmonic decomposition and a sampling scheme flag. 
     * The default resolution if the bandlimit of the input spherical harmonic transform.
     * @param flm the spherical harmonic decomposition
     * @param sc the flag of the sampling scheme
     * @return a new instance for either a MW or a Healpix map
     */
    static public PixelizedMap fromHarmonics(SphericalHarmonicTransform flm, SamplingScheme sc) {
        int defaultNside = flm.getBandlimit();
        return fromHarmonics(flm, sc, defaultNside);
    }
    /**
     * Factory method to create a Pixelized map from a wavelet decomposition, a sampling scheme flag and the resolution.
     * @param wavelets the axisymmetric wavelet decomposition
     * @param sc the flag of the sampling scheme
     * @param nside the Healpix resolution, only required if the map is Healpix
     * @return a new instance for either a MW or a Healpix map
     */
    static public PixelizedMap fromWavelets(AxisymmetricWaveletTransform wavelets, SamplingScheme sc, int nside) {
        if (sc == SamplingScheme.MW){
            if (wavelets.getReality()) {
                return new MWRealMap(wavelets);
            } else {
                return new MWCmplxMap(wavelets);
            }
        } else {
            return new HealpixRealMap(wavelets, nside);
        }
    };
    /**
     * Simplified factory method to create a Pixelized map from a wavelet decomposition and a sampling scheme flag. 
     * The default resolution if the bandlimit of the input wavelet transform.
     * @param wavelets the axisymmetric wavelet decomposition
     * @param sc the flag of the sampling scheme
     * @return a new instance for either a MW or a Healpix map
     */
     static public PixelizedMap fromWavelets(AxisymmetricWaveletTransform wavelets, SamplingScheme sc) {
        int defaultNside = wavelets.getBandlimit();
        return fromWavelets(wavelets, sc, defaultNside);
    }
     /**
      * Abstract method to compute the spherical harmonic decomposition of the map. 
      * The implementation depends on the map type.
      * @param flm the spherical harmonic decomposition
      */
    public abstract void computeHarmonics(Pointer<ComplexDouble> flm);
};

/**
 * Private class corresponding to a complex MW signal on the sphere.
 */
class MWCmplxMap extends PixelizedMap {
    /**
     * Private structure pointing to the data, i.e. the array of complex doubles on the sphere.
     */
    private Pointer<ComplexDouble> mapValues; 
    /**
     * Getter for the data structure.
     * @return the array of complex doubles
     */
    public Pointer<ComplexDouble> getMap() {
        return mapValues;
    }
    /**
     * Private constructor to construct a complex MW map from a spherical harmonic decomposition (using S2Let native methods).
     * @param flm the spherical harmonic decomposition
     */
    protected MWCmplxMap(SphericalHarmonicTransform flm) {
        super(SamplingScheme.MW, flm.getBandlimit(), false);
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1));
        s2let_mw_alm2map(this.mapValues, flm.getCoefficients(), resolution);
    }
    /**
     * Private constructor to construct a complex MW map from a wavelet decomposition (using S2Let native methods).
     * @param wavelets the axisymmetric wavelet decomposition
     */
    protected MWCmplxMap(AxisymmetricWaveletTransform wavelets) {
        super(SamplingScheme.MW, wavelets.getBandlimit(), false);
        this.mapValues = Pointer.allocateArray(ComplexDouble.class, resolution * (2*resolution-1));
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, resolution * resolution);
        wavelets.reconstructHarmonics(flm);
        s2let_mw_alm2map(this.mapValues, flm, resolution);
    }
    /**
     * Method to return the spherical harmonic coefficients of the map (using S2Let native methods).
     * @param flm a pointer to the structure you want to fill with the spherical harmonic coefficients
     */
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        s2let_mw_map2alm(flm, this.mapValues, resolution);
    };
}

/**
 * Private class corresponding to a real MW signal on the sphere.
 */
class MWRealMap extends PixelizedMap {
    /**
     * Private structure pointing to the data, i.e. the array of doubles on the sphere.
     */
    private Pointer<Double> mapValues; 
    /**
     * Getter for the data structure.
     * @return the array of doubles
     */
    public Pointer<Double> getMap() {
        return mapValues;
    }
    /**
     * Private constructor to construct a real MW map from a spherical harmonic decomposition (using S2Let native methods).
     * @param flm the spherical harmonic decomposition
     */
    protected MWRealMap(SphericalHarmonicTransform flm) {
        super(SamplingScheme.MW, flm.getBandlimit(), true);
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        s2let_mw_alm2map_real(this.mapValues, flm.getCoefficients(), resolution);
    }
    /**
     * Private constructor to construct a real MW map from a wavelet decomposition (using S2Let native methods).
     * @param wavelets the axisymmetric wavelet decomposition
     */
    protected MWRealMap(AxisymmetricWaveletTransform wavelets) {
        super(SamplingScheme.MW, wavelets.getBandlimit(), true);
        this.mapValues = Pointer.allocateArray(Double.class, resolution * (2*resolution-1));
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, resolution * resolution);
        wavelets.reconstructHarmonics(flm);
        s2let_mw_alm2map_real(this.mapValues, flm, resolution);
    }
    /**
     * Method to return the spherical harmonic coefficients of the map (using S2Let native methods).
     * @param flm a pointer to the structure you want to fill with the spherical harmonic coefficients
     */
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        s2let_mw_map2alm_real(flm, this.mapValues, resolution);
    };

}

/**
 * Private class corresponding to a real Healpix signal on the sphere.
 */
class HealpixRealMap extends PixelizedMap {
    /**
     * Private structure pointing to the data, i.e. the array of doubles on the sphere.
     */
    private Pointer<Double> mapValues; 
    /**
     * Getter for the data structure.
     * @return the array of doubles
     */
    public Pointer<Double> getMap() {
        return mapValues;
    }
    
    /**
     * Private constructor to construct a real Healpix map from a spherical harmonic decomposition (using S2Let native methods) at a given resolution.
     * @param flm the spherical harmonic decomposition
     * @param res the nside Healpix resolution
     */
    protected HealpixRealMap(SphericalHarmonicTransform flm, int res) {
        super(SamplingScheme.HEALPIX, res, true);
        this.mapValues = Pointer.allocateArray(Double.class, 12*res*res);
        s2let_hpx_alm2map_real(this.mapValues, flm.getCoefficients(), this.resolution, flm.getBandlimit());
    }
    /**
     * Private constructor to construct a real Healpix map from a wavelet decomposition (using S2Let native methods) at a given resolution.
     * @param wavelets the axisymmetric wavelet decomposition
     * @param res the nside Healpix resolution
     */
    protected HealpixRealMap(AxisymmetricWaveletTransform wavelets, int res) {
        super(SamplingScheme.HEALPIX, res, true);
        this.mapValues = Pointer.allocateArray(Double.class, 12*res*res);
        int bandlimit = wavelets.getBandlimit();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        wavelets.reconstructHarmonics(flm);
        s2let_hpx_alm2map_real(this.mapValues, flm, res, bandlimit);
    }
    /**
     * Method to return the spherical harmonic coefficients of the map (using S2Let native methods).
     * @param flm a pointer to the structure you want to fill with the spherical harmonic coefficients
     */
    public void computeHarmonics(Pointer<ComplexDouble> flm){
        int bandlimit = (int)(Math.sqrt(flm.getValidElements()));
        s2let_hpx_map2alm_real(flm, this.mapValues, resolution, bandlimit);
    };
}


