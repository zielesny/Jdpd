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
package de.gnwi.jdpd.samples.interactions.dpdCutoff1;

import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.parameters.InteractionDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.RandomAdderGroup;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Particle pair DPD random force calculator that uses 
 * ParticlePairInteractionCalculator to determine particle pair interactions 
 * that operate within a defined cut-off length of 1.0.
 * The passed particle x,y,z-position arrays are NOT changed. 
 * The following particle arrays are changed during accumulation process ("+-=" 
 * operations):
 * ParticleArrays.getF_x
 * ParticleArrays.getF_y
 * ParticleArrays.getF_z
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary.
 * 
 * @author Achim Zielesny
 */
public class ParticlePairDpdForceRandomCutoff1Calculator extends ParticlePairInteractionCalculator implements IParticlePairForceCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Numeric constants
     */
    private static final double ONE = 1.0;
    // </editor-fold>
    //
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
    public ParticlePairDpdForceRandomCutoff1Calculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
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
        this.simulationLogger.appendMethodCall("ParticlePairDpdForceRandomCutoff1Calculator.Constructor: FULL");
        // </editor-fold>
    }

    /**
     * Constructor that clones aParticlePairInteractionCalculator 
 (see ParticlePairInteractionCalculator for details)
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairDpdForceRandomCutoff1Calculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        super(aParticlePairInteractionCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairDpdForceRandomCutoff1Calculator.Constructor WITH aParticlePairInteractionCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected overridden methods">
    /**
     * Returns if force is assigned to f
     * 
     * @return True: Force is assigned to f, false: Otherwise
     */
    @Override
    public boolean isFAccumulation() {
        return true;
    }
    
    /**
     * Returns if force is assigned to ftwo
     * 
     * @return True: Force is assigned to ftwo, false: Otherwise
     */
    @Override
    public boolean isFtwoAccumulation() {
        return false;
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
        this.calculateParticlePairInteraction(
            aParticleIndex_i, 
            aParticleIndex_j, 
            aRij_x,
            aRij_y,
            aRij_z,
            aRij_Square,
            tmpRij,
            aRandomAdderGroup,
            aParameters
        );
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
        // Support quantities
        final ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        final InteractionDescription tmpInteractionDescription = aParameters.getInteractionDescription();
        final double tmpRandomValue;
        if (tmpInteractionDescription.isGaussianRandomDpdForce()) {
            tmpRandomValue = aRandomAdderGroup.getRandomNumberGenerator().nextGaussian();
        } else {
            tmpRandomValue = aRandomAdderGroup.getRandomNumberGenerator().nextZeroMeanUnitVarianceDouble();
        }
        final double tmpFactor = tmpInteractionDescription.getDpdSigmaDivRootTimeStepLength() * tmpRandomValue * (ONE/aRij - ONE);
        // Random DPD force
        final double tmpFij_x = tmpFactor * aRij_x;
        final double tmpFij_y = tmpFactor * aRij_y;
        final double tmpFij_z = tmpFactor * aRij_z;
        // Add to particle forces
        // NOTE: ParticlePairInteractionCalculator parallelization avoids collisions 
        //       of access of particle indices thus a lock is NOT necessary
        final double[] tmpF_x = tmpParticleArrays.getF_x();
        final double[] tmpF_y = tmpParticleArrays.getF_y();
        final double[] tmpF_z = tmpParticleArrays.getF_z();
        tmpF_x[aParticleIndex_i] += tmpFij_x;
        tmpF_x[aParticleIndex_j] -= tmpFij_x;
        tmpF_y[aParticleIndex_i] += tmpFij_y;
        tmpF_y[aParticleIndex_j] -= tmpFij_y;
        tmpF_z[aParticleIndex_i] += tmpFij_z;
        tmpF_z[aParticleIndex_j] -= tmpFij_z;
    }
    // </editor-fold>
    
}
