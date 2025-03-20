
package org.s2let;

import org.bridj.*;
import java.util.Random;

/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */

/**
 * Public class for storing a spherical harmonic transform, i.e. its harmonic coefficients and the relevant parameters.
 */
public class SphericalHarmonicTransform {
    /**
     * The spherical harmonic coefficients, an array of complex doubles.
     */
    private final Pointer<ComplexDouble> coefficients;
    /**
     * The bandlimit of the decomposition, i.e. the number of ell multipoles (size of flm = bandlimit^2).
     */
    private final int bandlimit;
    /** 
     * Reality flag: true if the complex coefficients in fact correcpond to a real signal on the sphere. 
     * This is true iff f(l,-m) = (-1)^m * CONG(f(l,m)).
     */
    private final boolean reality;
    /**
     * Getter for the data structure of the spherical harmonic coefficients.
     * @return the harmonic coefficients, i.e. an array of complex doubles
     */
    public Pointer<ComplexDouble> getCoefficients() {
        return coefficients;
    }
    /**
     * Getter for the bandlimit of the decomposition.
     * @return the bandlimit
     */
    public int getBandlimit() {
        return this.bandlimit;
    }
    /**
     * Getter for the reality flag.
     * @return true iff f(l,-m) = (-1)^m * CONG(f(l,m))
     */
    public boolean getReality() {
        return this.reality;
    }
    /**
     * Internal method to detect if the reality flag, i.e. if f(l,-m) = (-1)^m * CONG(f(l,m)) 
     * and thus if the decomposition corresponds to a real signal on the sphere.
     * @return true if f(l,-m) = (-1)^m * CONG(f(l,m)), false otherwise.
     */
    private boolean detectReality() {
        Pointer<Double> pd = coefficients.as(Double.class);
    	for (int el = 0; el < bandlimit; el++) {
            int em = 0;
            int i = el*el + el + em;
            double val = pd.getDoubleAtOffset(i * 2 * 8 + 8);
            if ( val != 0.0 ){
                return(false);
            }
            for (em = 1; em <= el; em++) {
                i = el*el + el + em;
                double val1 = pd.getDoubleAtOffset(i * 2 * 8);
                double val2 = pd.getDoubleAtOffset(i * 2 * 8 + 8);
                int iop = el*el + el - em;
                double val3 = pd.getDoubleAtOffset(iop * 2 * 8);
                double val4 = pd.getDoubleAtOffset(iop * 2 * 8 + 8);
                if ( val3 != Math.pow(-1.0,em) * val1 || val4 != - Math.pow(-1.0,em+1) * val2 ){
                    return(false);
                }
            }
        }
        return(true);
    }
    /**
     * Factory method to construct a spherical harmonic decomposition of given bandlimit filled with zeros.
     * @param bandlimit the bandlimit of the decomposition
     * @return a new instance of an empty spherical harmonic decomposition of given bandlimit
     */
    static public SphericalHarmonicTransform zeros(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        return new SphericalHarmonicTransform(flm);
    }
    /**
     * Factory method to construct a spherical harmonic decomposition of given bandlimit filled with random complex numbers.
     * @param bandlimit the bandlimit of the decomposition
     * @return a new instance of a random (complex) spherical harmonic decomposition of given bandlimit
     */
    static public SphericalHarmonicTransform random(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
	Random randomGenerator = new Random();
	Pointer<Double> pd = flm.as(Double.class);
    	for (int i = 0; i < bandlimit*bandlimit; i++) {
            pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
            pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
        }
        return new SphericalHarmonicTransform(flm);
    }
    /**
     * Factory method to construct a real spherical harmonic decomposition of given bandlimit filled with random complex numbers 
     * (but corresponding to a real signal, i.e. with f(l,-m) = (-1)^m * CONG(f(l,m))).
     * @param bandlimit the bandlimit of the decomposition
     * @return a new instance of a random (real) spherical harmonic decomposition of given bandlimit
     */
    static public SphericalHarmonicTransform randomReal(int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        Random randomGenerator = new Random();
	Pointer<Double> pd = flm.as(Double.class);
    	for (int el = 0; el < bandlimit; el++) {
            int em = 0;
            int i = el*el + el + em;
            pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
            pd.setDoubleAtOffset(i * 2 * 8 + 8, 0.0);
            for (em = 1; em <= el; em++) {
                i = el*el + el + em;
                pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
                pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
                int iop = el*el + el - em;
                pd.setDoubleAtOffset(iop * 2 * 8, Math.pow(-1.0,em) * pd.getDoubleAtOffset(i * 2 * 8));
                pd.setDoubleAtOffset(iop * 2 * 8 + 8, Math.pow(-1.0,em+1) * pd.getDoubleAtOffset(i * 2 * 8 + 8));
            }
        }
        return new SphericalHarmonicTransform(flm);
    }
    /**
     * Method for comparing the current decomposition with another spherical harmonic decomposition of identical band limit,
     * return the maximum absolute difference (comparing complex numbers two by two).
     * @param f2 another spherical harmonic decomposition of same bandlimit
     * @return the maximum absolute difference
     */
    public double maxAbsoluteDifferenceWith(Pointer<ComplexDouble> f2) {
        Pointer<ComplexDouble> f1 = coefficients;
        double val = 0.0;
        Pointer<Double> v1 = f1.as(Double.class);
        Pointer<Double> v2 = f2.as(Double.class);
        int len = Math.min( (int)f1.getValidElements(), (int)f2.getValidElements());
        for (int i = 0; i < len; i++) {
            val = Math.max( val, Math.max( v1.getDoubleAtOffset(i * 2 * 8) - v2.getDoubleAtOffset(i * 2 * 8), 
                    v1.getDoubleAtOffset(i * 2 * 8 + 8) - v2.getDoubleAtOffset(i * 2 * 8 + 8) ) );
        }
        return val;
    }
    /**
     * Method for comparing the current decomposition with another spherical harmonic decomposition of identical band limit,
     * return the maximum absolute difference (comparing complex numbers two by two).
     * @param f2 another spherical harmonic decomposition of same bandlimit
     * @return the maximum absolute difference
     */
    public double maxAbsoluteDifferenceWith(SphericalHarmonicTransform f2) {
        return maxAbsoluteDifferenceWith(f2.getCoefficients());
    }
    /**
     * Private constructor to create a new instance of spherical harmonic decomposition from the data (array of complex numbers) only.
     * @param flm the input data, i.e. Pointer to an array of complex doubles.
     */
    private SphericalHarmonicTransform(Pointer<ComplexDouble> flm) {
        this.coefficients = flm;
        this.bandlimit = (int)Math.sqrt(flm.getValidElements());
        this.reality = detectReality();
    }
    /**
     * Factory method to create an instance of spherical harmonic transform from the complex coefficients.
     * @param flm the spherical harmonic coefficients
     * @return a new instance of SphericalHarmonicTransform
     */
    static public SphericalHarmonicTransform fromHarmonics(Pointer<ComplexDouble> flm) {
        return new SphericalHarmonicTransform(flm);
    }
    /**
     * Factory method to create an instance of spherical harmonic transform from a pixelized map.
     * The spherical harmonic coefficients will be calculated in the map itself, using native methods from S2LET.
     * @param f a pixelized map, i.e. either a MW or a Healpix map
     * @return a new instance of SphericalHarmonicTransform
     */
    static public SphericalHarmonicTransform fromMap(PixelizedMap f) {
        int bandlimit = f.getResolution();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f.computeHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }  
    /**
     * Factory method to create an instance of spherical harmonic transform from a pixelized map at a given bandlimit.
     * The spherical harmonic coefficients will be calculated in the map itself, using native methods from S2LET.
     * @param f a pixelized map, i.e. either a MW or a Healpix map
     * @param bandlimit the bandlimit for the decomposition to create (important for Healpix maps)
     * @return a new instance of SphericalHarmonicTransform
     */ 
    static public SphericalHarmonicTransform fromMap(PixelizedMap f, int bandlimit) {
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f.computeHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }
    /**
     * Factory method to create an instance of spherical harmonic transform from a wavelet transform.
     * The spherical harmonic coefficients will be calculated in the map itself, using native methods from S2LET.
     * @param f_wav a axisymmetric wavelet transform
     * @return a new instance of SphericalHarmonicTransform
     */
    static public SphericalHarmonicTransform fromWavelets(AxisymmetricWaveletTransform f_wav) {     
        int bandlimit = f_wav.getBandlimit();
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);
        f_wav.reconstructHarmonics(flm);
        return new SphericalHarmonicTransform(flm);
    }
    
}
