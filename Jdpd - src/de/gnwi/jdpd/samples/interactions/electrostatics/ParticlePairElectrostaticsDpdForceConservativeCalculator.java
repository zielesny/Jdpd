/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpd.samples.interactions.electrostatics;

import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.Electrostatics;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.RandomAdderGroup;
import de.gnwi.jdpd.utilities.Utils;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.util.FastMath;

/**
 * Particle pair DPD electrostatics conservative force calculator
 * type that uses ParticlePairInteractionCalculator to determine particle pair 
 * interactions that operate within a defined cut-off length.
 * The passed particle x,y,z-position arrays are NOT changed.
 * The following particle arrays are changed during accumulation process ("+-=" 
 * operations):
 * ParticleArrays.getF_x
 * ParticleArrays.getF_y
 * ParticleArrays.getF_z
 * NOTE: ParticlePairInteractionCalculator parallelisation avoids collisions 
 * of access of particle indices thus a lock is NOT necessary
 * 
 * @author Achim Zielesny
 */
public class ParticlePairElectrostaticsDpdForceConservativeCalculator extends ParticlePairInteractionCalculator implements IParticlePairForceCalculator {

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
    public ParticlePairElectrostaticsDpdForceConservativeCalculator(
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
        this.simulationLogger.appendMethodCall("ParticlePairElectrostaticsDpdForceCalculator.Constructor: FULL");
        // </editor-fold>
    }

    /**
     * Constructor that clones aParticlePairInteractionCalculator (see 
     * ParticlePairInteractionCalculator for details)
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairElectrostaticsDpdForceConservativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        super(aParticlePairInteractionCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairElectrostaticsDpdForceCalculator.Constructor WITH aParticlePairInteractionCalculator");
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
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters,
        ParticlePairDistanceParameters aParticlePairDistanceParameters
    ) {
        // Support quantities
        final ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        final int[] tmpChargedParticleIndices = tmpParticleArrays.getChargedParticleIndices(); 
        final double[] tmpCharges = tmpParticleArrays.getCharges();
        // NOTE: ParticlePairElectrostaticsDpdForceConservativeCalculator is called with
        //       aParameters.getParticleArrays().getChargedParticles_r_x() etc.
        //       so that particle indices of method MUST be properly re-converted with 
        //       tmpParticleArrays.getChargedParticleIndices()
        final double tmpParticleCharge1 = tmpCharges[tmpChargedParticleIndices[aParticleIndex_i]];
        final double tmpParticleCharge2 = tmpCharges[tmpChargedParticleIndices[aParticleIndex_j]];
        if (tmpParticleCharge1 == 0.0 || tmpParticleCharge2 == 0.0) {
            throw new IllegalStateException("ParticlePairElectrostaticsDpdForceConservativeCalculator.calculateParticlePairInteraction: A charged particle does NOT have a charge.");
        }
        final Electrostatics tmpElectrostatics = aParameters.getInteractionDescription().getElectrostatics();
        final double tmpRij = FastMath.sqrt(aRij_Square);
        // Force = - dEpot/dRij
        double tmpForce = Utils.getElectrostaticsDpdForce(tmpParticleCharge1, tmpParticleCharge2, tmpElectrostatics, tmpRij);
        // Damping
        if (tmpRij <= tmpElectrostatics.getDampingDistance()) {
            tmpForce *= tmpElectrostatics.getDampingFactor();
        }
        // Maximum force
        if (tmpForce > tmpElectrostatics.getMaximumAbsoluteForceValue()) {
            tmpForce = tmpElectrostatics.getMaximumAbsoluteForceValue();
        } else if (tmpForce < -tmpElectrostatics.getMaximumAbsoluteForceValue()) {
            tmpForce = -tmpElectrostatics.getMaximumAbsoluteForceValue();
        }
        // Electrostatics force
        final double tmpFij_x = tmpForce * aRij_x / tmpRij;
        final double tmpFij_y = tmpForce * aRij_y / tmpRij;
        final double tmpFij_z = tmpForce * aRij_z / tmpRij;
        // Add to particle forces
        // NOTE: ParticlePairInteractionCalculator parallelization avoids collisions 
        //       of access of particle indices thus a lock is NOT necessary
        final double[] tmpF_x = tmpParticleArrays.getF_x();
        final double[] tmpF_y = tmpParticleArrays.getF_y();
        final double[] tmpF_z = tmpParticleArrays.getF_z();
        // NOTE: ParticlePairElectrostaticsDpdForceConservativeCalculator is called with
        //       aParameters.getParticleArrays().getChargedParticles_r_x() etc.
        //       so that particle indices of method MUST be properly re-converted with 
        //       tmpParticleArrays.getChargedParticleIndices()
        tmpF_x[tmpChargedParticleIndices[aParticleIndex_i]] += tmpFij_x;
        tmpF_x[tmpChargedParticleIndices[aParticleIndex_j]] -= tmpFij_x;
        tmpF_y[tmpChargedParticleIndices[aParticleIndex_i]] += tmpFij_y;
        tmpF_y[tmpChargedParticleIndices[aParticleIndex_j]] -= tmpFij_y;
        tmpF_z[tmpChargedParticleIndices[aParticleIndex_i]] += tmpFij_z;
        tmpF_z[tmpChargedParticleIndices[aParticleIndex_j]] -= tmpFij_z;
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
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        double aRij,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters
    ) {
        // Method not defined
        throw new IllegalStateException("ParticlePairElectrostaticsDpdForceConservativeCalculator.calculateParticlePairInteraction: Overload NOT supported.");
    }
    // </editor-fold>
    
}
