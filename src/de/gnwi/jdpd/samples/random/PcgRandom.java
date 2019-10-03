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

/**
 * PCG RNG implementation
 * NOTE: Implementation is NOT thread-safe.
 * 
 * @author Achim Zielesny
 */
public class PcgRandom implements IRandom {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    private final double SQRT_12 = Math.sqrt(12.0);
    private final double A_HALF = 0.5;
    
    /**
     * Random number generator
     */
    private final Pcg32 pcg32RandomNumberGenerator;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * NOTE: Implementation is NOT thread-safe.
     * 
     * @param aSeed Seed value (greater/equal 0)
     * @param aNumberOfWarmUpSteps Number of warm-up steps (greater/equal 0)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public PcgRandom(int aSeed, int aNumberOfWarmUpSteps) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSeed < 0) {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: aSeed < 0.");
        }
        if (aNumberOfWarmUpSteps < 0) {
            throw new IllegalArgumentException("ApacheCommonsRandom.Constructor: aNumberOfWarmUpSteps < 0.");
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Native seed">
        this.pcg32RandomNumberGenerator = new Pcg32((long) aSeed, (long) aSeed + 1L);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Warm-up">
        for (int i = 0; i < aNumberOfWarmUpSteps; i++) {
            this.pcg32RandomNumberGenerator.nextDouble();
        }
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
        return this.SQRT_12 * (this.pcg32RandomNumberGenerator.nextDouble() - this.A_HALF);
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
        return this.pcg32RandomNumberGenerator.nextDouble();
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
        return this.pcg32RandomNumberGenerator.nextGaussian();
    }
    // </editor-fold>
    
}
