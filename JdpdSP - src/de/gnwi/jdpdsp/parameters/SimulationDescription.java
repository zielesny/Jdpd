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
package de.gnwi.jdpdsp.parameters;

import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;

/**
 * Time step parameters
 * 
 * @author Achim Zielesny
 */
public class SimulationDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Number of time steps
     */
    private final int timeStepNumber;

    /**
     * Time step length
     */
    private final float timeStepLength;

    /**
     * timeStepLength/2.0
     */
    private final float timeStepLengthHalf;

    /**
     * (timeStepLength * timeStepLength)/2.0
     */
    private final float timeStepLengthSquareHalf;

    /**
     * Time step frequency for output
     */
    private final int timeStepFrequencyForOutput;
    
    /**
     * Number of initial potential energy minimization steps
     */
    private final int initialPotentialEnergyMinimizationStepNumber;
    
    /**
     * Type of initial potential energy minimization (true: All forces, false: DPD force only)
     */
    private final boolean isInitialPotentialEnergyMinimizationWithAllForces;
    
    /**
     * Flag for initial potential energy minimization step output.
     * True: Potential energy minimization step output is generated, 
     * false: Otherwise (NO output)
     */
    private final boolean isInitialPotentialEnergyMinimizationStepOutput;

    /**
     * Periodic boundaries
     */
    private final PeriodicBoundaries periodicBoundaries;

    /**
     * Flag for use of DPD unit masses:
     * True : DPD masses of all particles are set to 1
     * False: The DPD mass of the most lightweight particle (often water) is set to 1. 
     *        The masses of all other particles are set in relation to their 
     *        molar mass ratios to the most lightweight particle.
     */
    private final boolean isDpdUnitMass;

    /**
     * Number of initial velocity scaling steps
     */
    private final int numberOfInitialVelocityScalingSteps;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aTimeStepNumber Number of time steps
     * @param aTimeStepLength Time step length
     * @param aTimeStepFrequencyForOutput Time step frequency for output
     * @param anInitialPotentialEnergyMinimizationStepNumber Number of initial potential energy minimization steps
     * @param anIsInitialPotentialEnergyMinimizationWithAllForces Type of initial potential energy minimization (true: All forces, false: DPD force only)
     * @param anIsInitialPotentialEnergyMinimizationStepOutput Flag for initial potential energy minimization step output
     * @param aPeriodicBoundaries Periodic boundaries
     * @param anIsDpdUnitMass Flag for use of DPD unit masses
     * @param aNumberOfInitialVelocityScalingSteps Number of initial velocity scaling steps
     */
    public SimulationDescription(
        int aTimeStepNumber,
        float aTimeStepLength,
        int aTimeStepFrequencyForOutput,
        int anInitialPotentialEnergyMinimizationStepNumber,
        boolean anIsInitialPotentialEnergyMinimizationWithAllForces,
        boolean anIsInitialPotentialEnergyMinimizationStepOutput,
        PeriodicBoundaries aPeriodicBoundaries,
        boolean anIsDpdUnitMass,
        int aNumberOfInitialVelocityScalingSteps
    ) {
        this.timeStepNumber = aTimeStepNumber;
        this.timeStepLength = aTimeStepLength;
        // timeStepLength/2.0
        this.timeStepLengthHalf = this.timeStepLength/2.0f;
        // (timeStepLength * timeStepLength)/2.0
        this.timeStepLengthSquareHalf = (this.timeStepLength * this.timeStepLength)/2.0f;
        this.timeStepFrequencyForOutput = aTimeStepFrequencyForOutput;
        this.initialPotentialEnergyMinimizationStepNumber = anInitialPotentialEnergyMinimizationStepNumber;
        this.isInitialPotentialEnergyMinimizationWithAllForces = anIsInitialPotentialEnergyMinimizationWithAllForces;
        this.isInitialPotentialEnergyMinimizationStepOutput = anIsInitialPotentialEnergyMinimizationStepOutput;
        this.periodicBoundaries = aPeriodicBoundaries;
        this.isDpdUnitMass = anIsDpdUnitMass;
        this.numberOfInitialVelocityScalingSteps = aNumberOfInitialVelocityScalingSteps;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Number of time steps
     * 
     * @return Number of time steps
     */
    public int getTimeStepNumber() {
        return this.timeStepNumber;
    }

    /**
     * Time step length
     * 
     * @return Time step length
     */
    public float getTimeStepLength() {
        return this.timeStepLength;
    }

    /**
     * timeStepLength/2.0
     * 
     * @return timeStepLength/2.0
     */
    public float getTimeStepLengthHalf () {
        return this.timeStepLengthHalf;
    }

    /**
     * (timeStepLength * timeStepLength)/2.0
     * 
     * @return (timeStepLength * timeStepLength)/2.0
     */
    public float getTimeStepLengthSquareHalf () {
        return this.timeStepLengthSquareHalf;
    }

    /**
     * Time step frequency for output
     * 
     * @return Time step frequency for output
     */
    public int getTimeStepFrequencyForOutput() {
        return this.timeStepFrequencyForOutput;
    }

    /**
     * Number of initial potential energy minimization steps
     * 
     * @return Number of initial potential energy minimization steps
     */
    public int getInitialPotentialEnergyMinimizationStepNumber() {
        return this.initialPotentialEnergyMinimizationStepNumber;
    }

    /**
     * Type of initial potential energy minimization (true: All forces, false: DPD force only)
     * 
     * @return Type of initial potential energy minimization (true: All forces, false: DPD force only)
     */
    public boolean isInitialPotentialEnergyMinimizationWithAllForces() {
        return this.isInitialPotentialEnergyMinimizationWithAllForces;
    }

    /**
     * Flag for initial potential energy minimization step output.
     * 
     * @return True: Potential energy minimization step output is generated, 
     * false: Otherwise (NO output)
     */
    public boolean isInitialPotentialEnergyMinimizationStepOutput() {
        return this.isInitialPotentialEnergyMinimizationStepOutput;
    }
    
    /**
     * Periodic boundaries
     * 
     * @return Periodic boundaries
     */
    public PeriodicBoundaries getPeriodicBoundaries() {
        return this.periodicBoundaries;
    }

    /**
     * Flag for use of DPD unit masses:
     * True : DPD masses of all particles are set to 1
     * False: The DPD mass of the most lightweight particle (often water) is set to 1. 
     *        The masses of all other particles are set in relation to their 
     *        molar mass ratios to the most lightweight particle.
     * 
     * @return Flag for use of DPD unit masses
     */
    public boolean isDpdUnitMass() {
        return this.isDpdUnitMass;
    }
    
    /**
     * Number of initial velocity scaling steps
     * 
     * @return Number of initial velocity scaling steps
     */
    public int getNumberOfInitialVelocityScalingSteps() {
        return this.numberOfInitialVelocityScalingSteps;
    }
    // </editor-fold>
    
}
