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


import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondChunkArrays;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import java.util.LinkedList;

/**
 * Interface for harmonic bond property calculator
 * 
 * @author Achim Zielesny
 */
public interface IHarmonicBondPropertyCalculator extends ICalculator {

    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Calculates bond properties.
     * 
     * @param aBondChunkArraysList List of bond chunk arrays
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    boolean calculateBondProperties(
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        Parameters aParameters);
    
    /**
     * Box size
     * 
     * @return Box size
     */
    BoxSize getBoxSize();
    
    /**
     * Periodic boundaries
     * 
     * @return Periodic boundaries
     */
    PeriodicBoundaries getPeriodicBoundaries();

    /**
     * Parallelisation info
     * 
     * @return Parallelisation info
     */
    ParallelizationInfo getParallelizationInfo();
    
    /**
     * Simulation logger
     * 
     * @return Simulation logger
     */
    ILogger getSimulationLogger();
    // </editor-fold>
    
}
