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
package de.gnwi.jdpd.samples.harmonicBonds;

import de.gnwi.jdpd.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpd.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.AdderGroup;
import org.apache.commons.math3.util.FastMath;

/**
 * Harmonic bond conservative force calculator
 * The following particle arrays are changed during accumulation process ("+-=" 
 * operations):
 * ParticleArrays.getF_x
 * ParticleArrays.getF_y
 * ParticleArrays.getF_z
 * NOTE: HarmonicBondPropertyCalculator parallelisation avoids collisions of access of 
 * particle indices thus a lock is NOT necessary.
 * 
 * @author Achim Zielesny
 */
public class HarmonicBondForceConservativeCalculator extends HarmonicBondPropertyCalculator implements IHarmonicBondForceCalculator {

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
    public HarmonicBondForceConservativeCalculator(
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
        this.simulationLogger.appendMethodCall("BondForceCalculator.Constructor");
        // </editor-fold>
    }
    
    /**
     * Constructor that clones aHarmonicBondPropertyCalculator
 (see HarmonicBondPropertyCalculator for details)
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public HarmonicBondForceConservativeCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        super(aHarmonicBondPropertyCalculator);
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("BondForceCalculator.Constructor WITH aHarmonicBondPropertyCalculator");
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
     * Bond force calculation
     * NOTE: HarmonicBondPropertyCalculator parallelisation guarantees for 
     * aDoubleAdder being unique in current working thread thus NO thread-safe 
     * implementation of double adder is necessary.
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
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aBondLength,
        double aForceConstant,
        HarmonicBond.HarmonicBondBehaviour aHarmonicBondBehaviour,
        AdderGroup anAdderGroup,
        Parameters aParameters) {
        final ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        double tmpRij = FastMath.sqrt(aRij_x * aRij_x + aRij_y * aRij_y + aRij_z * aRij_z);
        // NOTE: Sign of tmpDeviation is important for correct force direction, see below!
        double tmpDeviation = tmpRij - aBondLength;
        boolean tmpIsCalculation = false;
        switch (aHarmonicBondBehaviour) {
            case DEFAULT:
                tmpIsCalculation = true;
                break;
            case ATTRACTIVE:
                if (tmpDeviation > 0.0) {
                    tmpIsCalculation = true;
                }        
                break;
            case REPULSIVE:
                if (tmpDeviation < 0.0) {
                    tmpIsCalculation = true;
                }        
                break;
        }
        if (tmpIsCalculation) {
            double tmpForce = aForceConstant * tmpDeviation;
            // Bond force
            final double tmpFij_x = tmpForce * aRij_x / tmpRij;
            final double tmpFij_y = tmpForce * aRij_y / tmpRij;
            final double tmpFij_z = tmpForce * aRij_z / tmpRij;
            // Add to particle forces
            // NOTE: HarmonicBondPropertyCalculator parallelization avoids collisions 
            //       of access of particle indices thus a lock is NOT necessary
            final double[] tmpF_x = tmpParticleArrays.getF_x();
            final double[] tmpF_y = tmpParticleArrays.getF_y();
            final double[] tmpF_z = tmpParticleArrays.getF_z();
            tmpF_x[aParticleIndex_i] -= tmpFij_x;
            tmpF_x[aParticleIndex_j] += tmpFij_x;
            tmpF_y[aParticleIndex_i] -= tmpFij_y;
            tmpF_y[aParticleIndex_j] += tmpFij_y;
            tmpF_z[aParticleIndex_i] -= tmpFij_z;
            tmpF_z[aParticleIndex_j] += tmpFij_z;
        }
    }
    // </editor-fold>
    
}
