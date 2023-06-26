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
package de.gnwi.jdpdsp.samples.interactions.dpdCutoff1;

import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionPnhlnCalculator;
import de.gnwi.jdpdsp.parameters.InteractionDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.parameters.ParticleArrays;
import de.gnwi.jdpdsp.parameters.SimulationDescription;
import de.gnwi.jdpdsp.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.utilities.Factory;
import de.gnwi.jdpdsp.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import de.gnwi.jdpdsp.utilities.RandomAdderGroup;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.util.FastMath;

/**
 * Particle pair velocity update calculator for PNHLN integration type 
 * that uses ParticlePairInteractionCalculator to determine particle pair 
 * interactions that operate within a defined cut-off length of 1.0.
 * The passed particle x,y,z-position arrays are NOT changed.
 * The following particle arrays are changed during accumulation process:
 * ParticleArrays.getV_x
 * ParticleArrays.getV_y
 * ParticleArrays.getV_z
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary.
 * 
 * @author Achim Zielesny
 */
public class ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator extends ParticlePairInteractionCalculator implements IParticlePairInteractionPnhlnCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Numeric constants
     */
    private static final float ONE = 1.0f;
    private static final float A_HALF = 0.5f;
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
    public ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
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
        this.simulationLogger.appendMethodCall("ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator.Constructor: FULL");
        // </editor-fold>
    }

    /**
     * Constructor that clones aParticlePairInteractionCalculator 
 (see ParticlePairInteractionCalculator for details)
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        super(aParticlePairInteractionCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator.Constructor WITH aParticlePairInteractionCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected overridden methods">
    /**
     * Returns accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     * (No checks are performed)
     * 
     * @return Accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     */
    @Override
    public float getAccumulatedPnhln_G_AddersSum() {
        return this.calculatorUtils.getAccumulatedPnhln_G_AddersSum();
    }

    /**
     * Particle pair interaction calculation
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or float adder 
     * is necessary.
     * (No checks are performed)
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
        float aRij_x,
        float aRij_y,
        float aRij_z,
        float aRij_Square,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters,
        ParticlePairDistanceParameters aParticlePairDistanceParameters
    ) {
        final float tmpRij = (float) FastMath.sqrt(aRij_Square);
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
     * NO thread-safe implementation of random number generator or float adder 
     * is necessary.
     * (No checks are performed)
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
        float aRij_x,
        float aRij_y,
        float aRij_z,
        float aRij_Square,
        float aRij,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters
    ) {
        // Support quantities
        final ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        final InteractionDescription tmpInteractionDescription = aParameters.getInteractionDescription();
        final SimulationDescription tmpSimulationDescription = aParameters.getSimulationDescription();
        final float[] tmpV_x = tmpParticleArrays.getV_x();
        final float[] tmpV_y = tmpParticleArrays.getV_y();
        final float[] tmpV_z = tmpParticleArrays.getV_z();

        final float tmpFactor1_i;
        final float tmpFactor1_j;
        final float tmpM_ij;
        if (aParameters.getSimulationDescription().isDpdUnitMass()) {
            tmpFactor1_i = ONE;
            tmpFactor1_j = ONE;
            tmpM_ij = A_HALF;
        } else {
            float[] tmpDpdMasses = tmpParticleArrays.getDpdMasses();
            float tmpMass_i = tmpDpdMasses[aParticleIndex_i];
            float tmpMass_j = tmpDpdMasses[aParticleIndex_j];
            tmpFactor1_i = ONE/tmpMass_i;
            tmpFactor1_j = ONE/tmpMass_j;
            tmpM_ij = tmpMass_i * tmpMass_j /(tmpMass_i + tmpMass_j);
        }

        final float tmpFactor2 = tmpM_ij / aRij_Square;

        // Difference of particle velocities: Vector v_ij
        final float tmpVij_x = tmpV_x[aParticleIndex_i] - tmpV_x[aParticleIndex_j];
        final float tmpVij_y = tmpV_y[aParticleIndex_i] - tmpV_y[aParticleIndex_j];
        final float tmpVij_z = tmpV_z[aParticleIndex_i] - tmpV_z[aParticleIndex_j];
        final float tmpFactor3 = tmpVij_x * aRij_x + tmpVij_y * aRij_y + tmpVij_z * aRij_z;
        
        final float tmpFactor4 = ONE - aRij - aRij + aRij_Square;
        final float tmpFactor5 = tmpSimulationDescription.getTimeStepLengthHalf() * tmpFactor2 * tmpFactor3 * ((float) FastMath.exp(-tmpInteractionDescription.getPnhln_Ksi() * tmpFactor4 * tmpSimulationDescription.getTimeStepLength() / tmpM_ij) - ONE);

        // Update v_i
        tmpV_x[aParticleIndex_i] += tmpFactor1_i * tmpFactor5 * aRij_x;
        tmpV_y[aParticleIndex_i] += tmpFactor1_i * tmpFactor5 * aRij_y;
        tmpV_z[aParticleIndex_i] += tmpFactor1_i * tmpFactor5 * aRij_z;
        // Update v_j
        tmpV_x[aParticleIndex_j] -= tmpFactor1_j * tmpFactor5 * aRij_x;
        tmpV_y[aParticleIndex_j] -= tmpFactor1_j * tmpFactor5 * aRij_y;
        tmpV_z[aParticleIndex_j] -= tmpFactor1_j * tmpFactor5 * aRij_z;
        
        // Calculate G if specified
        if (tmpInteractionDescription.isPnhln_G_Calculation()) {
            // New difference of particle velocities: Vector v_ij
            final float tmpVnew_ij_x = tmpV_x[aParticleIndex_i] - tmpV_x[aParticleIndex_j];
            final float tmpVnew_ij_y = tmpV_y[aParticleIndex_i] - tmpV_y[aParticleIndex_j];
            final float tmpVnew_ij_z = tmpV_z[aParticleIndex_i] - tmpV_z[aParticleIndex_j];
            final float tmpFactor6 = tmpVnew_ij_x * aRij_x + tmpVnew_ij_y * aRij_y + tmpVnew_ij_z * aRij_z;
            aRandomAdderGroup.getPnhln_G_Adder().add(tmpFactor4 * (ONE / aRij_Square * tmpFactor6 * tmpFactor6 - tmpInteractionDescription.getTemperature() / tmpM_ij));
        }
    }
    // </editor-fold>
    
}
