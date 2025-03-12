
package org.s2let;
import static org.s2let.bindings.S2letLibrary.*;
import org.bridj.*;


/**
 * @author      Boris Leistedt <boris.leistedt @ gmail.com>
 */

/**
 * Class to encapsulate the data and information of a axisymmetric wavelet transform.
 */
public class AxisymmetricWaveletTransform {
    /**
     * The content of the decomposition, i.e. spherical harmonic coefficients of the wavelet scales.
     */
    protected Pointer<ComplexDouble> f_wav_lm, f_scal_lm;
    /** 
     * The bandlimit of the decomposition, i.e. of the spherical harmonic coefficients of the wavelets.
     */
    protected final int bandlimit;
    /**
     * The wavelet parameter for the transform, i.e. the parameter for the tiling of harmonic space.
     */
    protected final int waveletParameter;
    /**
     * The first scale to use for the wavelets.
     */
    protected final int firstScale;
    /**
     * The number of wavelets, which depends on the bandlimit and the wavelet parameter.
     */
    protected final int scaleCount;
    /**
     * Reality flag, <code>true</code> if the spherical harmonic coefficients in fact all correspond to real signals (wavelets) on the sphere.
     */
    protected final boolean reality;
    /**
     * Multiresolution flag, <code>true</code> if the wavelets are calculated at different resolutions corresponding to their bandlimits. 
     */
    protected final boolean multiresolution; 
    /**
     * Getter for the multiresolution flag.
     * @return the multiresolution flag
     */
    public boolean getMultiresolution() {
        return multiresolution;
    }
    /**
     * Getter for the bandlimit
     * @return the bandlimit of the decomposition
     */
    public int getBandlimit() {
        return bandlimit;
    }
    /**
     * Getter for the first scale
     * @return the first scale of the decomposition
     */
    public int getFirstScale() {
        return firstScale;
    }
    /**
     * Getter for the wavelet parameter
     * @return the wavelet parameter of the decomposition
     */
    public int getWaveletParameter() {
        return waveletParameter;
    }
    /**
     * Getter for the scale count
     * @return the number of scales of the decomposition
     */
    public int getScaleCount() {
        return scaleCount;
    }
    /**
     * Getter for the reality flag
     * @return the reality flag of the decomposition
     */
    public boolean getReality() {
        return this.reality;
    }
    /**
     * Method for computing the bandlimit of the Nth scale of the decomposition (useful in the multiresolution case).
     * @param n the index of the scale
     * @return the bandlimit of the scale
     */
    public int getNthScaleBandlimit(int n) {
        return computeNthScaleBandlimit(n);
    }
    /**
     * Private for computing the bandlimit of the Nth scale of the decomposition (useful in the multiresolution case).
     * @param n the index of the scale
     * @return the bandlimit of the scale
     */
    private int computeNthScaleBandlimit(int n) {
        return (int)Math.ceil(Math.pow(waveletParameter, n+1));
    }
    /**
     * Method for computing the total size of the multiresolution wavelet transform.
     * @return the size of the transform in spherical harmonic space
     */
    private int getMultiresTotalSize_lm() {
        int tempBandlimit, total = 0;
        for(int j = firstScale; j <= scaleCount; j++){
            tempBandlimit = Math.min(getNthScaleBandlimit(j), bandlimit);
            total += tempBandlimit * tempBandlimit;
	}
        return total;
    }
    /**
     * Private constructor to create a new axisymmetric wavelet transform through an exact tiling of harmonic space,
     * from an input spherical harmonic transform (using S2LET native methods).
     * @param flm the input spherical harmonic transform
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param multires the multiresolution flag, i.e. <code>true</code> if all wavelet scales need to be constructed at the resolution corresponding to their bandlimit
     */
    private AxisymmetricWaveletTransform(SphericalHarmonicTransform flm, int waveletParameter, int firstScale, boolean multires) {
        this.reality = flm.getReality();
        this.bandlimit = flm.getBandlimit();
        this.scaleCount = s2let_j_max(bandlimit, waveletParameter);
        this.firstScale = firstScale;
        this.waveletParameter = waveletParameter;
        this.multiresolution = multires;
    	Pointer<Double> wav_lm = Pointer.allocateArray(Double.class, (scaleCount + 1) * bandlimit);
    	Pointer<Double> scal_lm = Pointer.allocateArray(Double.class, bandlimit);
    	s2let_axisym_wav_lm(wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        if(!multires){
            f_wav_lm = Pointer.allocateArray(ComplexDouble.class, (scaleCount + 1 - firstScale) * bandlimit * bandlimit);
            f_scal_lm = Pointer.allocateArray(ComplexDouble.class, bandlimit * bandlimit);   
            s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }else{
            f_wav_lm = Pointer.allocateArray(ComplexDouble.class, getMultiresTotalSize_lm());
            int bl = Math.min(computeNthScaleBandlimit(firstScale-1), bandlimit);
            f_scal_lm = Pointer.allocateArray(ComplexDouble.class, bl * bl);
            s2let_axisym_wav_analysis_multires_lm(f_wav_lm, f_scal_lm, flm.getCoefficients(), wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }
    }
    /**
     * Factory method to construct a axisymmetric wavelet transform from a pixelized map by doing the spherical harmonic transform and tiling the harmonic space to construct the wavelets.
     * @param f a pixelized map, either in MW or Healpix format
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param multires the multiresolution flag, i.e. <code>true</code> if all wavelet scales need to be constructed at the resolution corresponding to their bandlimit
     * @return a new instance of the spherical harmonic transform
     */
    static public AxisymmetricWaveletTransform fromMWMap(PixelizedMap f, int waveletParameter, int firstScale, boolean multires) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.fromMap(f);
        return new AxisymmetricWaveletTransform(flm, waveletParameter, firstScale, multires);
    }
    /**
     * Factory method to construct a axisymmetric wavelet transform from a pixelized map by doing the spherical harmonic transform and tiling the harmonic space to construct the wavelets.
     * @param f a pixelized map, either in MW or Healpix format
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param bandlimit the bandlimit of the spherical harmonic decomposition
     * @return a new instance of the spherical harmonic transform
     */
    static public AxisymmetricWaveletTransform fromHealpixMap(PixelizedMap f, int waveletParameter, int firstScale, int bandlimit) {
        SphericalHarmonicTransform flm = SphericalHarmonicTransform.fromMap(f, bandlimit);
        return new AxisymmetricWaveletTransform(flm, waveletParameter, firstScale, false);
    }
    /**
     * Factory method to create a new axisymmetric wavelet transform through an exact tiling of harmonic space,
     * from an input spherical harmonic transform (using S2LET native methods).
     * @param flm the input spherical harmonic transform
     * @param waveletParameter the wavelet parameter for the wavelet transform
     * @param firstScale the first scale to be used in the wavelet transform
     * @param multires the multiresolution flag, i.e. <code>true</code> if all wavelet scales need to be constructed at the resolution corresponding to their bandlimit
     * @return a new instance of the spherical harmonic transform
     */
    static public AxisymmetricWaveletTransform fromHarmonics(SphericalHarmonicTransform flm, int waveletParameter, int firstScale, boolean multires) {
        return new AxisymmetricWaveletTransform(flm, waveletParameter, firstScale, multires);
    }
    /**
     * Get the spherical harmonic transform the the Nth wavelet in the decomposition
     * @param scale the scale of interest
     * @return the spherical harmonic decomposition of the wavelet scale
     */
    public SphericalHarmonicTransform getNthScaleHarmonics(int scale) {
        int bl;
        if(multiresolution){
            bl = getNthScaleBandlimit(scale);
        } else {
            bl = this.bandlimit;
        }
        Pointer<ComplexDouble> flm = f_wav_lm.next(scale-firstScale).validElements(bl * bl);
        return SphericalHarmonicTransform.fromHarmonics(flm);
    }
    /**
     * Get the map (MW or Healpix) of the Nth wavelet in the decomposition
     * @param scale the scale of interest
     * @param sc the flag of the sampling scheme
     * @param nside resolution of the map
     * @return a pixelized map corresponding to the wavelet scale
     */
    public PixelizedMap computeNthScaleMap(int scale, SamplingScheme sc, int nside) {
        return PixelizedMap.fromHarmonics(getNthScaleHarmonics(scale), sc, nside);
    }
    /**
     * Get the map (MW or Healpix) of the Nth wavelet in the decomposition
     * @param scale the scale of interest
     * @param sc the flag of the sampling scheme
     * @param nside resolution of the map
     * @return a pixelized map corresponding to the wavelet scale
     */
    public PixelizedMap computeNthScaleMap(int scale, SamplingScheme sc) {
        return PixelizedMap.fromHarmonics(getNthScaleHarmonics(scale), sc);
    }
    /**
     * Reconstruct the spherical harmonics of the signal corresponding to the wavelet transform
     * @param flm the spherical harmonic transform
     */
    public void reconstructHarmonics(Pointer<ComplexDouble> flm){
        Pointer<Double> wav_lm = Pointer.allocateArray(Double.class, (scaleCount + 1) * bandlimit);
    	Pointer<Double> scal_lm = Pointer.allocateArray(Double.class, bandlimit);
    	s2let_axisym_wav_lm(wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        if(!multiresolution){
            s2let_axisym_wav_synthesis_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }else{
            s2let_axisym_wav_synthesis_multires_lm(flm, f_wav_lm, f_scal_lm, wav_lm, scal_lm, waveletParameter, bandlimit, firstScale);
        }
    }
}