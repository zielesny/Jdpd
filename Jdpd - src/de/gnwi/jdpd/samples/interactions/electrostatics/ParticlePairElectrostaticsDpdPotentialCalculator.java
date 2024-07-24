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
 * Particle pair DPD electrostatics potential calculator that uses 
 * ParticlePairInteractionCalculator to determine particle pair interactions 
 * that operate within a defined cut-off length.
 * NOTE: NO particle related array is changed.
 * 
 * @author Achim Zielesny
 */
public class ParticlePairElectrostaticsDpdPotentialCalculator extends ParticlePairInteractionCalculator {

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
    public ParticlePairElectrostaticsDpdPotentialCalculator(
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
        this.simulationLogger.appendMethodCall("ParticlePairElectrostaticsDpdPotentialCalculator.Constructor: FULL");
        // </editor-fold>
    }

    /**
     * Constructor that clones aParticlePairInteractionCalculator 
 (see ParticlePairInteractionCalculator for details)
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairElectrostaticsDpdPotentialCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        super(aParticlePairInteractionCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairElectrostaticsDpdPotentialCalculator.Constructor WITH aParticlePairInteractionCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected overridden methods">
    /**
     * Particle pair interaction calculation
     * Particle pair interaction calculation
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or double adder 
     * is necessary.
     * NOTE: aRandomAdderGroup is used for potential value accumulation.
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
            throw new IllegalStateException("ParticlePairElectrostaticsDpdPotentialCalculator.calculateParticlePairInteraction: A charged particle does NOT have a charge.");
        }
        final Electrostatics tmpElectrostatics = aParameters.getInteractionDescription().getElectrostatics();
        final double tmpRij = FastMath.sqrt(aRij_Square);
        // Force = - dEpot/dRij
        double tmpForce = Utils.getElectrostaticsDpdForce(tmpParticleCharge1, tmpParticleCharge2, tmpElectrostatics, tmpRij);
        // Damping
        if (tmpRij <= tmpElectrostatics.getDampingDistance()) {
            tmpForce *= tmpElectrostatics.getDampingFactor();
        }
        // Set potential energy according to |F| (for |F| = 0 there is no potential energy)
        if (tmpForce > tmpElectrostatics.getMaximumAbsoluteForceValue() || tmpForce < -tmpElectrostatics.getMaximumAbsoluteForceValue()) {
            // <editor-fold defaultstate="collapsed" desc="|F| > Fmax: Set Epot,max according to Fmax">
            double tmpPotentialEnergy;
            double tmpRijFmax;
            double tmpCoupledChargeProduct = tmpElectrostatics.getElectrostaticsCoupling() * tmpParticleCharge1 * tmpParticleCharge2;
            if (tmpElectrostatics.getEffectiveExponent() == 2.0) {
                tmpRijFmax = FastMath.sqrt(FastMath.abs(tmpCoupledChargeProduct)/tmpElectrostatics.getMaximumAbsoluteForceValue());
                tmpPotentialEnergy = tmpCoupledChargeProduct / tmpRijFmax;
            } else {
                tmpRijFmax = FastMath.pow(FastMath.abs(tmpCoupledChargeProduct)/tmpElectrostatics.getMaximumAbsoluteForceValue(), 1.0 / tmpElectrostatics.getEffectiveExponent());
                tmpPotentialEnergy = tmpCoupledChargeProduct / FastMath.pow(tmpRijFmax, tmpElectrostatics.getEffectiveExponent() - 1.0);
            }
            if (tmpRij <= tmpElectrostatics.getDampingDistance()) {
                tmpPotentialEnergy *= tmpElectrostatics.getDampingFactor();
            }
            aRandomAdderGroup.getPotentialEnergyAdder().add(tmpPotentialEnergy);
            // Pressure diagonal terms
            double tmpForceTerm = tmpElectrostatics.getMaximumAbsoluteForceValue() * FastMath.signum(tmpCoupledChargeProduct) / tmpRij;
            if (tmpRij <= tmpElectrostatics.getDampingDistance()) {
                tmpForceTerm *= tmpElectrostatics.getDampingFactor();
            }
            aRandomAdderGroup.getPressureXAdder().add(tmpForceTerm * aRij_x * aRij_x);
            aRandomAdderGroup.getPressureYAdder().add(tmpForceTerm * aRij_y * aRij_y);
            aRandomAdderGroup.getPressureZAdder().add(tmpForceTerm * aRij_z * aRij_z);
            // </editor-fold>
        } else if (tmpForce != 0.0) {
            // <editor-fold defaultstate="collapsed" desc="|F| <= Fmax and |F| != 0">
            double tmpPotentialEnergy = Utils.getElectrostaticsDpdPotentialEnergy(tmpParticleCharge1, tmpParticleCharge2, tmpElectrostatics, tmpRij);
            aRandomAdderGroup.getPotentialEnergyAdder().add(tmpPotentialEnergy);
            // Pressure diagonal terms
            aRandomAdderGroup.getPressureXAdder().add(tmpForce * aRij_x * aRij_x);
            aRandomAdderGroup.getPressureYAdder().add(tmpForce * aRij_y * aRij_y);
            aRandomAdderGroup.getPressureZAdder().add(tmpForce * aRij_z * aRij_z);
            // </editor-fold>
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
        throw new IllegalStateException("ParticlePairElectrostaticsDpdPotentialCalculator.calculateParticlePairInteraction: Overload NOT supported.");
    }
    // </editor-fold>
    
}
