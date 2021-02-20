/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2021  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.interfaces;

import de.gnwi.jdpd.accumulators.ParticleForceMagnitudeAccumulator;
import de.gnwi.jdpd.accumulators.PotentialAccumulator;

/**
 * Interface for time step calculator
 * 
 * @author Achim Zielesny
 */
public interface ITimeStepCalculator {
    
    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Calculates single time step
     * 
     * @param aCurrentTimeStep Current time step
     */
    void calculate(int aCurrentTimeStep);
    
    /**
     * Potential accumulator
     * 
     * @return Potential accumulator
     */
    PotentialAccumulator getPotentialAccumulator();

    /**
     * Particle force magnitude accumulator
     * 
     * @return Particle force magnitude accumulator
     */
    ParticleForceMagnitudeAccumulator getParticleForceMagnitudeAccumulator();
    
    /**
     * Executor services shutdown
     */
    void shutdownExecutorServices();
    // </editor-fold>
    
}
