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
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.RandomAdderGroup;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Particle pair DPD dissipative force calculator for SCMVV integration type 
 * that uses ParticlePairInteractionCalculator to determine particle pair 
 * interactions that operate within a defined cut-off length of 1.0.
 * The following global particle arrays are used (and NOT changed) in addition 
 * to the passed particle x,y,z-position arrays (which are also NOT changed):
 * ParticleArrays.getV_x
 * ParticleArrays.getV_y
 * ParticleArrays.getV_z
 * The following particle arrays are changed during accumulation process ("+-=" 
 * operations):
 * ParticleArrays.getFtwo_x
 * ParticleArrays.getFtwo_y
 * ParticleArrays.getFtwo_z
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary.
 * 
 * @author Achim Zielesny
 */
public class ParticlePairScmvvDpdForceDissipativeCutoff1Calculator extends ParticlePairInteractionCalculator implements IParticlePairForceCalculator {

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
    public ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(
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
    public ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
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
        return false;
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
        final double tmpRijInversMinusOne = ONE/aRij - 1.0;
        
        // Difference of particle velocities: Vector v_ij
        // IMPORTANT: Use tmpParticleArrays.getV_x() etc., NOT tmpParticleArrays.getVnew_x()
        final double[] tmpV_x = tmpParticleArrays.getV_x();
        final double[] tmpV_y = tmpParticleArrays.getV_y();
        final double[] tmpV_z = tmpParticleArrays.getV_z();
        final double tmpVij_x = tmpV_x[aParticleIndex_i] - tmpV_x[aParticleIndex_j];
        final double tmpVij_y = tmpV_y[aParticleIndex_i] - tmpV_y[aParticleIndex_j];
        final double tmpVij_z = tmpV_z[aParticleIndex_i] - tmpV_z[aParticleIndex_j];
        final double tmpFactor_v_dot_r = -aParameters.getInteractionDescription().getDpdGamma() * tmpRijInversMinusOne * tmpRijInversMinusOne * (tmpVij_x * aRij_x + tmpVij_y * aRij_y + tmpVij_z * aRij_z);
        
        // Dissipative DPD force
        final double tmpFtwo_ij_x = tmpFactor_v_dot_r * aRij_x;
        final double tmpFtwo_ij_y = tmpFactor_v_dot_r * aRij_y;
        final double tmpFtwo_ij_z = tmpFactor_v_dot_r * aRij_z;

        // IMPORTANT: Add to particle forces ftwo
        // NOTE: ParticlePairInteractionCalculator parallelization avoids collisions 
        //       of access of particle indices thus a lock is NOT necessary
        final double[] tmpFtwo_x = tmpParticleArrays.getFtwo_x();
        final double[] tmpFtwo_y = tmpParticleArrays.getFtwo_y();
        final double[] tmpFtwo_z = tmpParticleArrays.getFtwo_z();
        tmpFtwo_x[aParticleIndex_i] += tmpFtwo_ij_x;
        tmpFtwo_x[aParticleIndex_j] -= tmpFtwo_ij_x;
        tmpFtwo_y[aParticleIndex_i] += tmpFtwo_ij_y;
        tmpFtwo_y[aParticleIndex_j] -= tmpFtwo_ij_y;
        tmpFtwo_z[aParticleIndex_i] += tmpFtwo_ij_z;
        tmpFtwo_z[aParticleIndex_j] -= tmpFtwo_ij_z;
        
        // <editor-fold defaultstate="collapsed" desc="SCMVV logging">
        if (this.simulationLogger.isLogLevel(ILogger.SCMVV) && this.simulationLogger.isScmvvInformationAccumulation()) {
            this.simulationLogger.getScmvvFdissParticleIndexPairCounter().incrementAndGet();
            this.simulationLogger.getScmvvFdissRij_x_Adder().add(aRij_x);
            this.simulationLogger.getScmvvFdissAbsRij_x_Adder().add(Math.abs(aRij_x));
            this.simulationLogger.getScmvvFdissVij_x_Adder().add(tmpVij_x);
            this.simulationLogger.getScmvvFdissAbsVij_x_Adder().add(Math.abs(tmpVij_x));
            this.simulationLogger.getScmvvFdissRijVij_x_Adder().add(tmpVij_x * aRij_x);
            this.simulationLogger.getScmvvFdissAbsRijVij_x_Adder().add(Math.abs(tmpVij_x * aRij_x));
            this.simulationLogger.getScmvvFdissVdotR_Adder().add(tmpVij_x * aRij_x + tmpVij_y * aRij_y + tmpVij_z * aRij_z);
            this.simulationLogger.getScmvvFdissGammaFactor_Adder().add(aParameters.getInteractionDescription().getDpdGamma() * tmpRijInversMinusOne * tmpRijInversMinusOne);
        }
        // </editor-fold>
    }
    // </editor-fold>
    
}
