package org.s2let;

import static org.s2let.bindings.S2letLibrary.*;
import org.bridj.*;

import org.junit.*;
import static org.junit.Assert.*;

import static java.lang.Math.*;
import java.util.Random;

public class S2letLibraryTest {
        /**
         * Generate random complex harmonic coefficients.
         * @param flm
         * @param L 
         */
	public void random_flm(Pointer<ComplexDouble> flm, int L) {
            Random randomGenerator = new Random();
            Pointer<Double> pd = flm.as(Double.class);
            for (int i = 0; i < L*L; i++) {
                    pd.setDoubleAtOffset(i * 2 * 8, randomGenerator.nextDouble());
                    pd.setDoubleAtOffset(i * 2 * 8 + 8, randomGenerator.nextDouble());
            }
	}
        /**
         * Generate random complex harmonic coefficients corresponding to a real signal.
         * @param flm
         * @param L 
         */
        public void random_flm_real(Pointer<ComplexDouble> flm, int L) {
            Random randomGenerator = new Random();
            Pointer<Double> pd = flm.as(Double.class);
            for (int el = 0; el < L; el++) {
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
	}

	// Test native function s2let_j_max
	@Test
	public void test_s2let_j_max() {
		test_s2let_j_max(64, 2);
		test_s2let_j_max(128, 3);
		test_s2let_j_max(256, 4);
	}
	public void test_s2let_j_max(int L, int B) {
		assertEquals((int)ceil(log(L) / log(B)), s2let_j_max(L, B));
    }

    // Test native functions in harmonic space (full resolution axisym wavelet transform)
    @Test
    public void test_s2let_axisym_wav_lm() {
    	test_s2let_axisym_wav_lm(2, 4, 0);
    	test_s2let_axisym_wav_lm(2, 64, 1);
    	test_s2let_axisym_wav_lm(3, 128, 2);
    }
    public void test_s2let_axisym_wav_lm(int B, int L, int J_min) {
    	
    	Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, L * L);
    	random_flm(flm, L);

    	int J = s2let_j_max(L, B);
    	Pointer<Double> wav_lm = Pointer.allocateArray(Double.class, (J + 1) * L);
    	Pointer<Double> scal_lm = Pointer.allocateArray(Double.class, L);
    	s2let_axisym_wav_lm(wav_lm, scal_lm, B, L, J_min);

    	Pointer<ComplexDouble> f_wav_lm = Pointer.allocateArray(ComplexDouble.class, (J + 1 - J_min) * L * L);
    	Pointer<ComplexDouble> f_scal_lm = Pointer.allocateArray(ComplexDouble.class, L * L);
    	s2let_axisym_wav_analysis_lm(f_wav_lm, f_scal_lm, flm, wav_lm, scal_lm, B, L, J_min);

    	Pointer<ComplexDouble> flm_rec = Pointer.allocateArray(ComplexDouble.class, L * L);
    	s2let_axisym_wav_synthesis_lm(flm_rec, f_wav_lm, f_scal_lm, wav_lm, scal_lm, B, L, J_min);

    	Pointer<Double> v1 = flm.as(Double.class);
    	Pointer<Double> v2 = flm_rec.as(Double.class);
    	for (int i = 0; i < L*L; i++) {

    		assertEquals(v1.getDoubleAtOffset(i * 2 * 8), v2.getDoubleAtOffset(i * 2 * 8), 1e-8);
    		assertEquals(v1.getDoubleAtOffset(i * 2 * 8 + 8), v2.getDoubleAtOffset(i * 2 * 8 + 8), 1e-8);
    	}
    }

    @Test
    public void test_ssht_real() {
        test_ssht_real(4);
    }
    public void test_ssht_real(int L) {
        
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, L * L);
        Pointer<Double> f = Pointer.allocateArray(Double.class, L * (2*L-1));
        Pointer<ComplexDouble> flm_rec = Pointer.allocateArray(ComplexDouble.class, L * L);
       
        random_flm_real(flm, L);

        s2let_mw_alm2map_real(f, flm, L);
        s2let_mw_map2alm_real(flm_rec, f, L);

        Pointer<Double> v1 = flm.as(Double.class);
        Pointer<Double> v2 = flm_rec.as(Double.class);
        for (int i = 0; i < L*L; i++) {
            assertEquals(v1.getDoubleAtOffset(i * 2 * 8), v2.getDoubleAtOffset(i * 2 * 8), 1e-8);
            assertEquals(v1.getDoubleAtOffset(i * 2 * 8 + 8), v2.getDoubleAtOffset(i * 2 * 8 + 8), 1e-8);
        }
    }
    
    @Test
    public void test_ssht() {
        test_ssht(4);
        test_ssht(64);
        test_ssht(128);
    }
    public void test_ssht(int L) {
        
        Pointer<ComplexDouble> flm = Pointer.allocateArray(ComplexDouble.class, L * L);
        Pointer<ComplexDouble> f = Pointer.allocateArray(ComplexDouble.class, L * (2*L-1));
        Pointer<ComplexDouble> flm_rec = Pointer.allocateArray(ComplexDouble.class, L * L);
       
        random_flm(flm, L);

        s2let_mw_alm2map(f, flm, L);
        s2let_mw_map2alm(flm_rec, f, L);

        Pointer<Double> v1 = flm.as(Double.class);
        Pointer<Double> v2 = flm_rec.as(Double.class);
        for (int i = 0; i < L*L; i++) {
            assertEquals(v1.getDoubleAtOffset(i * 2 * 8), v2.getDoubleAtOffset(i * 2 * 8), 1e-8);
            assertEquals(v1.getDoubleAtOffset(i * 2 * 8 + 8), v2.getDoubleAtOffset(i * 2 * 8 + 8), 1e-8);
        }
    }

}