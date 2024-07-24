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
package de.gnwi.jdpdsp.parameters;

/**
 * Global parameters for simulation
 * 
 * @author Achim Zielesny
 */
public final class Parameters {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Restart info
     */
    private final RestartInfo restartInfo;

    /**
     * Particle arrays
     */
    private final ParticleArrays particleArrays;

    /**
     * Particle types
     */
    private final ParticleTypes particleTypes;    

    /**
     * Chemical system description
     */
    private final ChemicalSystemDescription chemicalSystemDescription;
    
    /**
     * Interaction description
     */
    private final InteractionDescription interactionDescription;
    
    /**
     * Simulation description
     */
    private final SimulationDescription simulationDescription;
    
    /**
     * Simulation counts
     */
    private final SimulationCounts simulationCounts;
    
    /**
     * Test objects
     */
    private final TestObjects testObjects;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aRestartInfo Restart info
     * @param aParticleArrays Particle arrays
     * @param aParticleTypes Particle types
     * @param aChemicalSystemDescription Chemical system description
     * @param anInteractionDescription Interaction description
     * @param aSimulationDescription Simulation description
     * @param aSimulationCounts Simulation counts
     * @param aTestObjects Test objects
     */
    public Parameters(
        RestartInfo aRestartInfo,
        ParticleArrays aParticleArrays,
        ParticleTypes aParticleTypes,    
        ChemicalSystemDescription aChemicalSystemDescription,
        InteractionDescription anInteractionDescription,
        SimulationDescription aSimulationDescription,
        SimulationCounts aSimulationCounts,
        TestObjects aTestObjects) {
        this.restartInfo = aRestartInfo;
        this.particleArrays = aParticleArrays;
        this.particleTypes = aParticleTypes;
        this.chemicalSystemDescription = aChemicalSystemDescription;
        this.interactionDescription = anInteractionDescription;
        this.simulationDescription = aSimulationDescription;
        this.simulationCounts = aSimulationCounts;
        this.testObjects = aTestObjects;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Returns if restart info is available
     * @return True: Restart info is available, false: Otherwise
     */
    public boolean hasRestartInfo() {
        return this.restartInfo != null;
    }
    
    /**
     * Returns maximum time step with velocity scaling
     * 
     * @return Maximum time step with velocity scaling
     */
    public int getMaximumTimeStepWithVelocityScaling() {
        if (this.hasRestartInfo()) {
            if (this.restartInfo.isVelocityInitialization()) {
                return this.restartInfo.getLastTimeStep() + this.simulationDescription.getNumberOfInitialVelocityScalingSteps();
            } else {
                return this.simulationDescription.getNumberOfInitialVelocityScalingSteps();
            }
        } else {
            return this.simulationDescription.getNumberOfInitialVelocityScalingSteps();
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Restart info
     * 
     * @return Restart info
     */
    public RestartInfo getRestartInfo() {
        return this.restartInfo;
    }

    /**
     * Particle arrays
     * 
     * @return Particle arrays
     */
    public ParticleArrays getParticleArrays() {
        return this.particleArrays;
    }

    /**
     * Particle types
     * 
     * @return Particle types
     */
    public ParticleTypes getParticleTypes() {
        return this.particleTypes;
    }    

    /**
     * Chemical system description
     * 
     * @return Chemical system description
     */
    public ChemicalSystemDescription getChemicalSystemDescription() {
        return this.chemicalSystemDescription;
    }
    
    /**
     * Interaction description
     * 
     * @return Interaction description
     */
    public InteractionDescription getInteractionDescription() {
        return this.interactionDescription;
    }
    
    /**
     * Simulation description
     * 
     * @return Simulation description
     */
    public SimulationDescription getSimulationDescription() {
        return this.simulationDescription;
    }
    
    /**
     * Simulation counts
     * 
     * @return Simulation counts
     */
    public SimulationCounts getSimulationCounts() {
        return this.simulationCounts;
    }
    
    /**
     * Test objects
     * 
     * @return Test objects
     */
    public TestObjects getTestObjects() {
        return this.testObjects;
    }
    // </editor-fold>

}
