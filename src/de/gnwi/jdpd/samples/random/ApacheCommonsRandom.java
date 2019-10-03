/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2019  Achim Zielesny (achim.zielesny@googlemail.com)
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
import org.apache.commons.rng.UniformRandomProvider;
import org.apache.commons.rng.sampling.distribution.ZigguratNormalizedGaussianSampler;
import org.apache.commons.rng.sampling.distribution.ContinuousSampler;
import org.apache.commons.rng.simple.RandomSource;

/**
 * Apache Commons RNG implementation
 * NOTE: Implementation is NOT thread-safe.
 * 
 * @author Achim Zielesny
 */
public class ApacheCommonsRandom implements IRandom {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    private final double SQRT_12 = Math.sqrt(12.0);
    private final double A_HALF = 0.5;
    
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
    // <editor-fold defaultstate="collapsed" desc="Constructor">
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
        // <editor-fold defaultstate="collapsed" desc="Native seed">
        boolean tmpIsIntArray = false;
        boolean tmpIsLongArray = false;
        int tmpArrayLength = -1;
        switch(aRandomSource) {
            case ISAAC:
                tmpIsIntArray = true;
                tmpArrayLength = 256;
                break;
            case JDK:
                tmpIsLongArray = true;
                tmpArrayLength = 1;
                break;
            case KISS:
                tmpIsIntArray = true;
                tmpArrayLength = 4;
                break;
            case MT:
                tmpIsIntArray = true;
                tmpArrayLength = 624;
                break;
            case MT_64:
                tmpIsLongArray = true;
                tmpArrayLength = 312;
                break;
            case MWC_256:
                tmpIsIntArray = true;
                tmpArrayLength = 257;
                break;
            case SPLIT_MIX_64:
                tmpIsLongArray = true;
                tmpArrayLength = 1;
                break;
            case TWO_CMRES:
                tmpIsIntArray = true;
                tmpArrayLength = 1;
                break;
            case WELL_1024_A:
                tmpIsIntArray = true;
                tmpArrayLength = 32;
                break;
            case WELL_19937_A:
                tmpIsIntArray = true;
                tmpArrayLength = 624;
                break;
            case WELL_19937_C:
                tmpIsIntArray = true;
                tmpArrayLength = 624;
                break;
            case WELL_44497_A:
                tmpIsIntArray = true;
                tmpArrayLength = 1391;
                break;
            case WELL_44497_B:
                tmpIsIntArray = true;
                tmpArrayLength = 1391;
                break;
            case WELL_512_A:
                tmpIsIntArray = true;
                tmpArrayLength = 16;
                break;
            case XOR_SHIFT_1024_S:
                tmpIsLongArray = true;
                tmpArrayLength = 16;
                break;
            default:
                throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: Unknown random source.");
        }
        if (tmpIsIntArray) {
            if (tmpArrayLength == 1) {
                this.randomNumberGenerator = RandomSource.create(aRandomSource, aSeed);
            } else {
                // Seed generation according to Takuji Nishimura and Makoto 
                // Matsumoto for MT19937 (improved initialization 2002/1/26)
                int[] tmpNativeIntegerSeedArray = new int[tmpArrayLength];
                tmpNativeIntegerSeedArray[0]= aSeed & 0xFFFFFFFF;
                for (int i = 1; i < tmpArrayLength; i++) {
                    tmpNativeIntegerSeedArray[i] = (1812433253 * (tmpNativeIntegerSeedArray[i - 1] ^ (tmpNativeIntegerSeedArray[i - 1] >> 30)) + i);
                    tmpNativeIntegerSeedArray[i] &= 0xFFFFFFFF;
                }
                this.randomNumberGenerator = RandomSource.create(aRandomSource, tmpNativeIntegerSeedArray);
            }
        } else if (tmpIsLongArray) {
            if (tmpArrayLength == 1) {
                this.randomNumberGenerator = RandomSource.create(aRandomSource, (long) aSeed);
            } else {
                // Seed generation according to Takuji Nishimura and Makoto 
                // Matsumoto for MT19937 (improved initialization 2002/1/26)
                long[] tmpNativeLongSeedArray = new long[tmpArrayLength];
                tmpNativeLongSeedArray[0] = ((long) aSeed) & 0xFFFFFFFFL;
                for (int i = 1; i < tmpArrayLength; i++) {
                    tmpNativeLongSeedArray[i] = (1812433253L * (tmpNativeLongSeedArray[i - 1] ^ (tmpNativeLongSeedArray[i - 1] >> 30)) + ((long) i));
                    tmpNativeLongSeedArray[i] &= 0xFFFFFFFFL;
                }            
                this.randomNumberGenerator = RandomSource.create(aRandomSource, tmpNativeLongSeedArray);
            }
        } else {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: Unknown random source.");
        }
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
    // </editor-fold>
    
}
