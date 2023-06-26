/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2023  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpdsp.interfaces;

/**
 * Interface for random number generator
 * 
 * @author Achim Zielesny
 */
public interface IRandom {

    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Uniform continuous pseudorandom with zero mean and unit variance
     * 
     * @return Uniform continuous pseudorandom with zero mean and unit variance
     */
    float nextZeroMeanUnitVarianceFloat();
    
    /**
     * Returns the next pseudorandom, uniformly distributed float value 
     * between 0.0 and 1.0 from this random number generator's sequence.
     * 
     * @return The next pseudorandom, uniformly distributed float value 
     * between 0.0 and 1.0 from this random number generator's sequence.
     */
    float nextFloat();
    
    /**
     * Returns the next pseudorandom, Gaussian ("normally") distributed float 
     * value with mean 0.0 and standard deviation 1.0 from this random number 
     * generator's sequence.
     * 
     * @return The next pseudorandom, Gaussian ("normally") distributed float 
     * value with mean 0.0 and standard deviation 1.0 from this random number 
     * generator's sequence.
     */
    float nextGaussian();
    // </editor-fold>
    
}
