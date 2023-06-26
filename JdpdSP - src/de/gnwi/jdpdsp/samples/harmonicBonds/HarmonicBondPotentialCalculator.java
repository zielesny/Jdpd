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
package de.gnwi.jdpdsp.samples.harmonicBonds;

import de.gnwi.jdpdsp.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import de.gnwi.jdpdsp.utilities.AdderGroup;
import org.apache.commons.math3.util.FastMath;

/**
 * Harmonic bond potential calculator
 * NOTE: NO particle related array is changed.
 * 
 * @author Achim Zielesny
 */
public class HarmonicBondPotentialCalculator extends HarmonicBondPropertyCalculator {

    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public HarmonicBondPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        super(
            aSimulationLogger, 
            aBoxSize, 
            aPeriodicBoundaries, 
            aParallelizationInfo
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("BondPotentialCalculator.Constructor");
        // </editor-fold>
    }
    
    /**
     * Constructor that clones aHarmonicBondPropertyCalculator
 (see HarmonicBondPropertyCalculator for details)
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public HarmonicBondPotentialCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        super(aHarmonicBondPropertyCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("BondPotentialCalculator.Constructor WITH aHarmonicBondPropertyCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected overridden methods">
    /**
     * Bond force calculation
     * NOTE: HarmonicBondPropertyCalculator parallelisation guarantees for 
     * aDoubleAdder being unique in current working thread thus NO thread-safe 
     * implementation of float adder is necessary.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aBondLength Bond length
     * @param aForceConstant Force constant
     * @param aHarmonicBondBehaviour Behaviour of harmonic bond
     * @param anAdderGroup Adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     */
    @Override
    protected void calculateBondProperty(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        float aRij_x,
        float aRij_y,
        float aRij_z,
        float aBondLength,
        float aForceConstant,
        HarmonicBond.HarmonicBondBehaviour aHarmonicBondBehaviour,
        AdderGroup anAdderGroup,
        Parameters aParameters
    ) {
        final float tmpRij_x_Square = aRij_x * aRij_x;
        final float tmpRij_y_Square = aRij_y * aRij_y;
        final float tmpRij_z_Square = aRij_z * aRij_z;
        final float tmpRij = (float) FastMath.sqrt(tmpRij_x_Square + tmpRij_y_Square + tmpRij_z_Square);
        final float tmpDeviation = tmpRij - aBondLength;
        boolean tmpIsCalculation = false;
        switch (aHarmonicBondBehaviour) {
            case DEFAULT:
                tmpIsCalculation = true;
                break;
            case ATTRACTIVE:
                if (tmpDeviation > 0.0f) {
                    tmpIsCalculation = true;
                }        
                break;
            case REPULSIVE:
                if (tmpDeviation < 0.0f) {
                    tmpIsCalculation = true;
                }        
                break;
        }
        if (tmpIsCalculation) {
            final float tmpTerm = aForceConstant * tmpDeviation;
            // Potential energy
            // NOTE: Factor 0.5 of harmonic spring potential is assumed to be already included in force constant k, i.e. k = 0.5 * k(real)
            anAdderGroup.getPotentialEnergyAdder().add(tmpTerm * tmpDeviation);
            // Pressure diagonal terms
            final float tmpForceTerm = - tmpTerm / tmpRij;
            anAdderGroup.getPressureXAdder().add(tmpForceTerm * tmpRij_x_Square);
            anAdderGroup.getPressureYAdder().add(tmpForceTerm * tmpRij_y_Square);
            anAdderGroup.getPressureZAdder().add(tmpForceTerm * tmpRij_z_Square);
        }
    }
    // </editor-fold>
        
}