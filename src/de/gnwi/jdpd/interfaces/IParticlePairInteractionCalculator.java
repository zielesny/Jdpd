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

import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Interface for particle pair interaction calculator
 * 
 * @author Achim Zielesny
 */
public interface IParticlePairInteractionCalculator extends ICalculator {

    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Calculates particle pair interactions with a distance smaller than 
     * cut-off length based on calculation mode.
     * NOTE: No checks are performed.
     * NOTE: This method scales linearly with the number of particles.
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters
     * @param aCalculationMode Calculation mode
     * @return True: Operation successful, false: Otherwise
     * @throws IllegalArgumentException Thrown if calculation mode is unknown
     */
    public boolean calculateParticlePairInteractions(
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aCalculationMode
    );

    /**
     * Sets activity of cache for ParticlePairDistanceParameterss
     * 
     * @param aValue True: Cache for ParticlePairDistanceParameterss is active, false: Otherwise
     */
    void setParticlePairDistanceParametersCacheActivity(boolean aValue);

    /**
     * Activity of cache for ParticlePairDistanceParameterss
     * 
     * @return True: Cache for ParticlePairDistanceParameterss is active, false: Otherwise
     */
    boolean getParticlePairDistanceParametersCacheActivity();
    
    /**
     * Sets cache for ParticlePairDistanceParameterss
     * 
     * @param aParticlePairDistanceParametersCache Cache for ParticlePairDistanceParameterss
     */
    void setParticlePairDistanceParametersCache(ParticlePairDistanceParameters[][] aParticlePairDistanceParametersCache);
    
    /**
     * Sets cache for ParticlePairDistanceParameterss
     * 
     * @param aParticlePairInteractionCalculator Particle pair interaction 
     * calculator whose cache is transferred
     */
    void setParticlePairDistanceParametersCache(IParticlePairInteractionCalculator aParticlePairInteractionCalculator);

    /**
     * Removes cache for ParticlePairDistanceParameterss
     */
    void removeParticlePairDistanceParametersCache();
    
    /**
     * Returns cache for ParticlePairDistanceParameterss
     * 
     * @return Cache for ParticlePairDistanceParameterss
     */
    ParticlePairDistanceParameters[][] getParticlePairDistanceParametersCache();
    
    /**
     * Set particle assignments (e.g. for call of method
 this.calculateCellBasedParticlePairInteractionsWithCurrentParticleCellAssignments())
 NOTE: NO checks are performed.
     * 
     * @param aHeadParticleIndexOfCellArray Head particle index of cell array
     * (retrieved with corresponding getter)
     * @param aNextParticleIndexOfCellArray Next particle index of cell array
     * (retrieved with corresponding getter)
     */
    public void setParticleCellAssignments(int[] aHeadParticleIndexOfCellArray, int[] aNextParticleIndexOfCellArray);

    /**
     * Set particle assignments (e.g. for call of method
     * this.calculateCellBasedParticlePairInteractionsWithCurrentParticleCellAssignments())
     * NOTE: NO checks are performed.
     * 
     * @param aParticlePairInteractionCalculator Particle pair interaction calculator whose 
     * head particle index of cell array and next particle index of cell array are 
     * transferred
     */
    public void setParticleCellAssignments(IParticlePairInteractionCalculator aParticlePairInteractionCalculator);

    /**
     * Head particle index of cell array
     * 
     * @return Head particle index of cell array or null if isCellBasedCalculation() 
     * is false
     */
    public int[] getHeadParticleIndexOfCellArray();

    /**
     * Next particle index of cell array
     * 
     * @return Next particle index of cell array or null if isCellBasedCalculation() 
     * is false
     */
    public int[] getNextParticleIndexOfCellArray();
            
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
     * Cut-off length for partitioning of the box
     * 
     * @return Cut-off length for partitioning of the box
     */
    double getCutOffLength();

    /**
     * Parallelisation info
     * 
     * @return Parallelisation info
     */
    ParallelizationInfo getParallelizationInfo();
    
    /**
     * Cell neighbours
     * 
     * @return Cell neighbours or null if isCellBasedCalculation() is false
     */
    int[][] getCellNeigbours();

    /**
     * Cell chunks that are parallelisation-safe
     * 
     * @return Cell chunks that are parallelisation-safe
     */
    int[][] getParallelizationSafeCellChunks();
    
    /**
     * Random number seed
     * 
     * @return Random number seed
     */
    AtomicInteger getRandomNumberSeed();
    
    /**
     * Simulation logger
     * 
     * @return Simulation logger
     */
    ILogger getSimulationLogger();
    
    /**
     * Factory for new objects
     * 
     * @return Factory for new objects
     */
    Factory getFactory();
    // </editor-fold>
    
}
