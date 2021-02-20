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
package de.gnwi.jdpd.samples.interactions.nearestNeighbor;

import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.RandomAdderGroup;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Particle pair nearest-neighbor calculator that uses 
 * ParticlePairInteractionCalculator to determine particle pair interactions 
 * that operate within a defined cut-off length of 1.0.
  * The following particle arrays are changed:
 * ParticleArrays.getNearestNeighborDistances
 * ParticleArrays.getNearestNeighborParticleIndices
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary.
* 
 * @author Achim Zielesny
 */
public class ParticlePairNearestNeighborCalculator extends ParticlePairInteractionCalculator {

    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aFactory Factory for new objects
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairNearestNeighborCalculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException 
    {
        super(
            aFactory,
            aSimulationLogger, 
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength,
            aParallelizationInfo, 
            aRandomNumberSeed
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairNearestNeighborCalculator.Constructor: FULL");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected overridden methods">
    /**
     * Particle pair interaction calculation
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or double adder 
     * is necessary.
     * NOTE: aRandomAdderGroup is used for potential value accumulation.
     * NOTE: No checks are performed.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aRij_Square Squared distance between particle i and j
     * @param aRandomAdderGroup Random adder group (NOTE: NOT thread-safe, see above)
     * @param aParameters Parameters (may be null)
     * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance ((may be null)
     */
    @Override
    protected void calculateParticlePairInteraction(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters,
        ParticlePairDistanceParameters aParticlePairDistanceParameters
    ) {
        final double tmpRij = Math.sqrt(aRij_Square);
        // Check if distance is small enough
        if (tmpRij <= aParameters.getChemicalSystemDescription().getNearestNeighborDistance()) {
            final boolean[] tmpIsNearestNeighborBaseParticleDeterminations = aParameters.getParticleArrays().isNearestNeighborBaseParticleDeterminations();
            // Check if nearest-neighbor determination for at least one of the particles
            if (tmpIsNearestNeighborBaseParticleDeterminations[aParticleIndex_i] || tmpIsNearestNeighborBaseParticleDeterminations[aParticleIndex_j]) {
                final int[] tmpMoleculeIndices = aParameters.getParticleArrays().getMoleculeIndices();
                // Check if both particles belong to DIFFERENT molecules
                if (tmpMoleculeIndices[aParticleIndex_i] != tmpMoleculeIndices[aParticleIndex_j]) {
                    if (tmpIsNearestNeighborBaseParticleDeterminations[aParticleIndex_i]) {
                        this.setArrays(
                            aParticleIndex_i,
                            aParticleIndex_j,
                            tmpRij,
                            aParameters.getParticleArrays().getNearestNeighborDistances(),
                            aParameters.getParticleArrays().getNearestNeighborParticleIndices()
                        );
                    }
                    if (tmpIsNearestNeighborBaseParticleDeterminations[aParticleIndex_j]) {
                        this.setArrays(
                            aParticleIndex_j,
                            aParticleIndex_i,
                            tmpRij,
                            aParameters.getParticleArrays().getNearestNeighborDistances(),
                            aParameters.getParticleArrays().getNearestNeighborParticleIndices()
                        );
                    }
                }
            }
        }
        // Add to particle pair distance parameters cache if necessary
        if (this.isParticlePairDistanceParametersCacheActive) {
            aParticlePairDistanceParameters.add(
                aParticleIndex_i, 
                aParticleIndex_j, 
                aRij_x,
                aRij_y,
                aRij_z,
                aRij_Square,
                tmpRij
            );
        }
    }

    /**
     * Particle pair interaction calculation
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or double adder 
     * is necessary.
     * NOTE: No checks are performed.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aRij_Square Squared distance between particle i and j
     * @param aRij Distance between particle i and j
     * @param aRandomAdderGroup Random adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     * @throws IllegalStateException Thrown if method is illegally called
     */
    @Override
    protected void calculateParticlePairInteraction(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        double aRij,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters
    ) {
        // Method not defined
        throw new IllegalStateException("ParticlePairNearestNeighborCalculator.calculateParticlePairInteraction: Overload NOT supported.");
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Sets arrays
     * NOTE: NO checks are performed
     * 
     * @param aSimulationLogger Simulation logger
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij Distance between particle i and j
     * @param aNearestNeighborDistances Distances to nearest (non-bonded) 
     * neighbor particles (of another molecule)
     * @param aNearestNeighborParticleIndices Indices of nearest (non-bonded) 
     * neighbor particles (of another molecule)
     */
    private void setArrays(
        int aParticleIndex_i,
        int aParticleIndex_j,
        double aRij,
        double[] aNearestNeighborDistances,
        int[] aNearestNeighborParticleIndices
    ) {
        if (aRij < aNearestNeighborDistances[aParticleIndex_i]) {
            aNearestNeighborDistances[aParticleIndex_i] = aRij;
            aNearestNeighborParticleIndices[aParticleIndex_i] = aParticleIndex_j;
        }
    }
    // </editor-fold>
    
}
