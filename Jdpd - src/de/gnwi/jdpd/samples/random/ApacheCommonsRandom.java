/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
 * 
 * Source code is available at <https://github.com/zielesny/Jdpd>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package de.gnwi.jdpd.samples.random;

import de.gnwi.jdpd.interfaces.IRandom;
import org.apache.commons.rng.JumpableUniformRandomProvider;
import org.apache.commons.rng.LongJumpableUniformRandomProvider;
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.simple.RandomSource;
import org.apache.commons.math3.util.FastMath;

/**
 * Apache Commons RNG implementation
 * NOTE: Implementation is NOT thread-safe.
 * 
 * @author Achim Zielesny
 */
public class ApacheCommonsRandom implements IRandom {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    private final double SQRT_12 = FastMath.sqrt(12.0);
    private final double A_HALF = 0.5;
    
    /**
     * Random source this RNG is constructed with
     */
    private final RandomSource randomSource;
    
    /**
     * Random number generator
     */
    private final UniformRandomProvider randomNumberGenerator;
    
    /**
     * Gaussian sampler
     */
    private final ContinuousSampler gaussianSampler;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * NOTE: Implementation is NOT thread-safe.
     * 
     * @param aRandomSource Random source
     * @param aSeed Seed value (greater/equal 0)
     * @param aNumberOfWarmUpSteps Number of warm-up steps (greater/equal 0)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ApacheCommonsRandom(RandomSource aRandomSource, int aSeed, int aNumberOfWarmUpSteps) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSeed < 0) {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: aSeed < 0.");
        }
        if (aNumberOfWarmUpSteps < 0) {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: aNumberOfWarmUpSteps < 0.");
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Initialize">
        this.randomSource = aRandomSource;
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Native seed">
        this.randomNumberGenerator = this.getRandomNumberGenerator(this.randomSource, aSeed);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Warm-up">
        for (int i = 0; i < aNumberOfWarmUpSteps; i++) {
            this.randomNumberGenerator.nextDouble();
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Gaussian sampler">
        // Old code Apache Commons RNG 1.0:
        // this.gaussianSampler = new BoxMullerGaussianSampler(this.randomNumberGenerator, 0.0, 1.0);
        this.gaussianSampler = new ZigguratNormalizedGaussianSampler(this.randomNumberGenerator);
        // </editor-fold>
    }

    /**
     * Constructor
     * NOTE: Implementation is NOT thread-safe.
     * 
     * @param aRandomSource Random source
     * @param aSeed Seed value (greater/equal 0)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ApacheCommonsRandom(RandomSource aRandomSource, int aSeed) {
        this(aRandomSource, aSeed, 0);
    }
    
    /**
     * Constructor
     * NOTE: Implementation is NOT thread-safe.
     * 
     * @param aRandomSource Random source
     * @param aRandomNumberGenerator Jumped random number generator
     */
    private ApacheCommonsRandom(RandomSource aRandomSource, UniformRandomProvider aRandomNumberGenerator) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aRandomNumberGenerator == null) {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: aRandomNumberGenerator is null.");
        }
        // </editor-fold>
        this.randomSource = aRandomSource;
        this.randomNumberGenerator = aRandomNumberGenerator;
        this.gaussianSampler = new ZigguratNormalizedGaussianSampler(this.randomNumberGenerator);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Uniform continuous pseudorandom with zero mean and unit variance
     * 
     * @return Uniform continuous pseudorandom with zero mean and unit variance
     */
    @Override
    public double nextZeroMeanUnitVarianceDouble() {
        return this.SQRT_12 * (this.randomNumberGenerator.nextDouble() - this.A_HALF);
    }
    
    /**
     * Returns the next pseudorandom, uniformly distributed double value 
     * between 0.0 and 1.0 from this random number generator's sequence.
     * 
     * @return The next pseudorandom, uniformly distributed double value 
     * between 0.0 and 1.0 from this random number generator's sequence.
     */
    @Override
    public double nextDouble() {
        return this.randomNumberGenerator.nextDouble();
    }
    
    /**
     * Returns the next pseudorandom, Gaussian ("normally") distributed double 
     * value with mean 0.0 and standard deviation 1.0 from this random number 
     * generator's sequence.
     * 
     * @return The next pseudorandom, Gaussian ("normally") distributed double 
     * value with mean 0.0 and standard deviation 1.0 from this random number 
     * generator's sequence.
     */
    @Override
    public double nextGaussian() {
        return this.gaussianSampler.sample();
    }
    
    /**
     * Returns if RNG is jumpable
     * @return True: RNG is jumpable, false: Otherwise
     */
    public boolean isJumpable() {
        return this.randomSource.isJumpable() || this.randomSource.isLongJumpable();
    }
    
    /**
     * Returns jumped RNG
     * 
     * @return Jumped RNG or null if RNG cannot jump
     */
    public ApacheCommonsRandom getJumpedRng() {
        if (this.randomSource.isLongJumpable()) {
            return new ApacheCommonsRandom(this.randomSource, ((LongJumpableUniformRandomProvider) this.randomNumberGenerator).longJump());
        } else if (this.randomSource.isJumpable()) {
            return new ApacheCommonsRandom(this.randomSource, ((JumpableUniformRandomProvider) this.randomNumberGenerator).jump());
        } else {
            return null;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Returns native integer seed array.
     * Seed generation according to Takuji Nishimura and Makoto Matsumoto 
     * for MT19937 (improved initialization 2002/1/26): init_genrand()
     * 
     * @param aSeed Seed value (greater/equal 0)
     * @param tmpArrayLength Array length
     * @return Native integer seed array
     */
    private int[] getNativeIntegerSeedArray(int aSeed, int tmpArrayLength) {
        int[] tmpNativeIntegerSeedArray = new int[tmpArrayLength];
        tmpNativeIntegerSeedArray[0]= aSeed & 0xFFFFFFFF;
        for (int i = 1; i < tmpArrayLength; i++) {
            tmpNativeIntegerSeedArray[i] = (1812433253 * (tmpNativeIntegerSeedArray[i - 1] ^ (tmpNativeIntegerSeedArray[i - 1] >> 30)) + i);
            tmpNativeIntegerSeedArray[i] &= 0xFFFFFFFF;
        }
        return tmpNativeIntegerSeedArray;
    }

    /**
     * Returns native long seed array.
     * Seed generation according to Takuji Nishimura and Makoto Matsumoto 
     * for MT19937-64 (2004/9/29 version): init_genrand64()
     * 
     * @param aSeed Seed value (greater/equal 0)
     * @param tmpArrayLength Array length
     * @return Native long seed array
     */
    private long[] getNativeLongSeedArray(int aSeed, int tmpArrayLength) {
        long[] tmpNativeLongSeedArray = new long[tmpArrayLength];
        tmpNativeLongSeedArray[0] = (long) aSeed;
        for (int i = 1; i < tmpArrayLength; i++) {
            tmpNativeLongSeedArray[i] = (6364136223846793005L * (tmpNativeLongSeedArray[i - 1] ^ (tmpNativeLongSeedArray[i - 1] >> 62)) + ((long) i));
        }            
        return tmpNativeLongSeedArray;
    }

    /**
     * Returns random number generator
     * 
     * @param aRandomSource Random source
     * @param aSeed Seed value (greater/equal 0)
     * @return Random number generator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    private UniformRandomProvider getRandomNumberGenerator (RandomSource aRandomSource, int aSeed) {
        int tmpArrayLength;
        UniformRandomProvider tmpRandomNumberGenerator;
        switch(aRandomSource) {
            case ISAAC:
                tmpArrayLength = 256;
                tmpRandomNumberGenerator = RandomSource.ISAAC.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case JDK:
                tmpRandomNumberGenerator = RandomSource.JDK.create((long) aSeed);
                break;
            case JSF_32: // new in ACRNG 1.4
                tmpRandomNumberGenerator = RandomSource.JSF_32.create(aSeed);
                break;
            case JSF_64: // new in ACRNG 1.4
                tmpRandomNumberGenerator = RandomSource.JSF_64.create((long) aSeed);
                break;
            case KISS:
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.KISS.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case MT:
                tmpArrayLength = 624;
                tmpRandomNumberGenerator = RandomSource.MT.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case MT_64:
                tmpArrayLength = 312;
                tmpRandomNumberGenerator = RandomSource.MT_64.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case MSWS: // new in ACRNG 1.4
                tmpArrayLength = 3;
                tmpRandomNumberGenerator = RandomSource.MSWS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case MWC_256:
                tmpArrayLength = 257;
                tmpRandomNumberGenerator = RandomSource.MWC_256.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case PCG_XSH_RR_32: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.PCG_XSH_RR_32.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case PCG_XSH_RS_32: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.PCG_XSH_RS_32.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case PCG_RXS_M_XS_64: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.PCG_RXS_M_XS_64.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case PCG_MCG_XSH_RR_32: // new in ACRNG 1.4
                tmpRandomNumberGenerator = RandomSource.PCG_MCG_XSH_RR_32.create((long) aSeed);
                break;
            case PCG_MCG_XSH_RS_32: // new in ACRNG 1.4
                tmpRandomNumberGenerator = RandomSource.PCG_MCG_XSH_RS_32.create((long) aSeed);
                break;
            case SFC_32: // new in ACRNG 1.4
                tmpArrayLength = 3;
                tmpRandomNumberGenerator = RandomSource.SFC_32.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case SFC_64: // new in ACRNG 1.4
                tmpArrayLength = 3;
                tmpRandomNumberGenerator = RandomSource.SFC_64.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case SPLIT_MIX_64:
                tmpRandomNumberGenerator = RandomSource.SPLIT_MIX_64.create((long) aSeed);
                break;
            case TWO_CMRES:
                tmpRandomNumberGenerator = RandomSource.TWO_CMRES.create(aSeed);
                break;
            case WELL_1024_A:
                tmpArrayLength = 32;
                tmpRandomNumberGenerator = RandomSource.WELL_1024_A.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case WELL_19937_A:
                tmpArrayLength = 624;
                tmpRandomNumberGenerator = RandomSource.WELL_19937_A.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case WELL_19937_C:
                tmpArrayLength = 624;
                tmpRandomNumberGenerator = RandomSource.WELL_19937_C.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case WELL_44497_A:
                tmpArrayLength = 1391;
                tmpRandomNumberGenerator = RandomSource.WELL_44497_A.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case WELL_44497_B:
                tmpArrayLength = 1391;
                tmpRandomNumberGenerator = RandomSource.WELL_44497_B.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case WELL_512_A:
                tmpArrayLength = 16;
                tmpRandomNumberGenerator = RandomSource.WELL_512_A.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XOR_SHIFT_1024_S_PHI: // new in ACRNG 1.4
                tmpArrayLength = 16;
                tmpRandomNumberGenerator = RandomSource.XOR_SHIFT_1024_S_PHI.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_64_S: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_64_S.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_64_SS: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_64_SS.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_128_PLUS: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_128_PLUS.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_128_SS: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_128_SS.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_128_PLUS: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_128_PLUS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_128_SS: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_128_SS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_256_PLUS: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_256_PLUS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_256_SS: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_256_SS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_512_PLUS: // new in ACRNG 1.4
                tmpArrayLength = 8;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_512_PLUS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_512_SS: // new in ACRNG 1.4
                tmpArrayLength = 8;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_512_SS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_128_PP: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_128_PP.create(this.getNativeIntegerSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_128_PP: // new in ACRNG 1.4
                tmpArrayLength = 2;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_128_PP.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_256_PP: // new in ACRNG 1.4
                tmpArrayLength = 4;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_256_PP.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_SHI_RO_512_PP: // new in ACRNG 1.4
                tmpArrayLength = 8;
                tmpRandomNumberGenerator = RandomSource.XO_SHI_RO_512_PP.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_1024_PP: // new in ACRNG 1.4
                tmpArrayLength = 16;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_1024_PP.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_1024_S: // new in ACRNG 1.4
                tmpArrayLength = 16;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_1024_S.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            case XO_RO_SHI_RO_1024_SS: // new in ACRNG 1.4
                tmpArrayLength = 16;
                tmpRandomNumberGenerator = RandomSource.XO_RO_SHI_RO_1024_SS.create(this.getNativeLongSeedArray(aSeed, tmpArrayLength));
                break;
            default:
                throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: Unknown random source.");
        }
        return tmpRandomNumberGenerator;
    }
    // </editor-fold>
    
}
