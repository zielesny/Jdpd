/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpdsp.samples.interactions.dpdCutoff1;

import de.gnwi.jdpdsp.utilities.Factory;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.parameters.ParticleArrays;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.parameters.InteractionDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import de.gnwi.jdpdsp.utilities.RandomAdderGroup;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.util.FastMath;

/**
 * Particle pair DPD full (i.e. conservative + random + dissipative) force 
 * calculator for SCMVV integration type that uses 
 * ParticlePairInteractionCalculator to determine particle pair interactions 
 * that operate within a defined cut-off length of 1.0.
 * The following global particle arrays are used (and NOT changed) in addition 
 * to the passed particle x,y,z-position arrays (which are also NOT changed):
 * ParticleArrays.getV_x
 * ParticleArrays.getV_y
 * ParticleArrays.getV_z
 * The following particle arrays are changed during accumulation process ("+-=" 
 * operations):
 * Conservative + Random part:
 * ParticleArrays.getF_x
 * ParticleArrays.getF_y
 * ParticleArrays.getF_z
 * Dissipative part:
 * ParticleArrays.getFtwo_x
 * ParticleArrays.getFtwo_y
 * ParticleArrays.getFtwo_z
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary.
 * 
 * @author Achim Zielesny
 */
public class ParticlePairScmvvDpdForceFullCutoff1Calculator extends ParticlePairInteractionCalculator implements IParticlePairForceCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Numeric constants
     */
    private static final float ONE = 1.0f;
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
    public ParticlePairScmvvDpdForceFullCutoff1Calculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        float aCutOffLength,
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
        this.simulationLogger.appendMethodCall("ParticlePairDpdRcut1ForceCalculator.Constructor: FULL");
        // </editor-fold>
    }

    /**
     * Constructor that clones aParticlePairInteractionCalculator 
 (see ParticlePairInteractionCalculator for details)
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairScmvvDpdForceFullCutoff1Calculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        super(aParticlePairInteractionCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairDpdRcut1ForceCalculator.Constructor WITH aParticlePairInteractionCalculator");
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
        return true;
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
        // Support quantities
        final ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        // <editor-fold defaultstate="collapsed" desc="Conservative + random force to f">
        final int[] tmpParticleTypeIndices = tmpParticleArrays.getParticleTypeIndices();
        final InteractionDescription tmpInteractionDescription = aParameters.getInteractionDescription();
        final float tmpRandomValue;
        if (tmpInteractionDescription.isGaussianRandomDpdForce()) {
            tmpRandomValue = aRandomAdderGroup.getRandomNumberGenerator().nextGaussian();
        } else {
            tmpRandomValue = aRandomAdderGroup.getRandomNumberGenerator().nextZeroMeanUnitVarianceFloat();
        }
        final float tmpRij = (float) FastMath.sqrt(aRij_Square);
        final float tmpRijInversMinusOne = ONE/tmpRij - 1.0f;
        final float tmpFactor1 = (tmpInteractionDescription.getAij()[tmpParticleTypeIndices[aParticleIndex_i]][tmpParticleTypeIndices[aParticleIndex_j]] + tmpInteractionDescription.getDpdSigmaDivRootTimeStepLength() * tmpRandomValue) * tmpRijInversMinusOne;
        // Conservative + Random DPD force
        final float tmpFij_x = tmpFactor1 * aRij_x;
        final float tmpFij_y = tmpFactor1 * aRij_y;
        final float tmpFij_z = tmpFactor1 * aRij_z;
        // IMPORTANT: Add to old particle forces f
        // NOTE: ParticlePairInteractionCalculator parallelization avoids collisions 
        //       of access of particle indices thus a lock is NOT necessary
        final float[] tmpF_x = tmpParticleArrays.getF_x();
        final float[] tmpF_y = tmpParticleArrays.getF_y();
        final float[] tmpF_z = tmpParticleArrays.getF_z();
        tmpF_x[aParticleIndex_i] += tmpFij_x;
        tmpF_x[aParticleIndex_j] -= tmpFij_x;
        tmpF_y[aParticleIndex_i] += tmpFij_y;
        tmpF_y[aParticleIndex_j] -= tmpFij_y;
        tmpF_z[aParticleIndex_i] += tmpFij_z;
        tmpF_z[aParticleIndex_j] -= tmpFij_z;
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Dissipative force to ftwo">
        // Difference of particle velocities: Vector v_ij
        // IMPORTANT: Use tmpParticleArrays.getV_x() etc., NOT tmpParticleArrays.getVnew_x()
        final float[] tmpV_x = tmpParticleArrays.getV_x();
        final float[] tmpV_y = tmpParticleArrays.getV_y();
        final float[] tmpV_z = tmpParticleArrays.getV_z();
        final float tmpVij_x = tmpV_x[aParticleIndex_i] - tmpV_x[aParticleIndex_j];
        final float tmpVij_y = tmpV_y[aParticleIndex_i] - tmpV_y[aParticleIndex_j];
        final float tmpVij_z = tmpV_z[aParticleIndex_i] - tmpV_z[aParticleIndex_j];
        final float tmpFactor_v_dot_r = -aParameters.getInteractionDescription().getDpdGamma() * tmpRijInversMinusOne * tmpRijInversMinusOne * (tmpVij_x * aRij_x + tmpVij_y * aRij_y + tmpVij_z * aRij_z);
        // Dissipative DPD force
        final float tmpFtwo_ij_x = tmpFactor_v_dot_r * aRij_x;
        final float tmpFtwo_ij_y = tmpFactor_v_dot_r * aRij_y;
        final float tmpFtwo_ij_z = tmpFactor_v_dot_r * aRij_z;
        // IMPORTANT: Add to particle forces ftwo
        // NOTE: ParticlePairInteractionCalculator parallelization avoids collisions 
        //       of access of particle indices thus a lock is NOT necessary
        final float[] tmpFtwo_x = tmpParticleArrays.getFtwo_x();
        final float[] tmpFtwo_y = tmpParticleArrays.getFtwo_y();
        final float[] tmpFtwo_z = tmpParticleArrays.getFtwo_z();
        tmpFtwo_x[aParticleIndex_i] += tmpFtwo_ij_x;
        tmpFtwo_x[aParticleIndex_j] -= tmpFtwo_ij_x;
        tmpFtwo_y[aParticleIndex_i] += tmpFtwo_ij_y;
        tmpFtwo_y[aParticleIndex_j] -= tmpFtwo_ij_y;
        tmpFtwo_z[aParticleIndex_i] += tmpFtwo_ij_z;
        tmpFtwo_z[aParticleIndex_j] -= tmpFtwo_ij_z;
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Add to particle pair distance parameters cache if necessary">
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
        // </editor-fold>
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
        // Method not defined
        throw new IllegalStateException("ParticlePairScmvvDpdForceFullCutoff1Calculator.calculateParticlePairInteraction: Overload NOT supported.");
    }
    // </editor-fold>
    
}
