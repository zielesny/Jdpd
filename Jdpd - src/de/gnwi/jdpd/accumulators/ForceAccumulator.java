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
package de.gnwi.jdpd.accumulators;

import de.gnwi.jdpd.samples.harmonicBonds.HarmonicBondChunkArrays;
import de.gnwi.jdpd.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * Force accumulator that uses (and changes) the particle related arrays
 * ParticleArrays.getF_x/getFtwo_x, 
 * ParticleArrays.getF_y/getFtwo_y and
 * ParticleArrays.getF_z/getFtwo_z
 * for accumulative ("+=" and "-=") operations.
 * 
 * @author Achim Zielesny
 */
public class ForceAccumulator {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Simulation logger
     */
    private final ILogger simulationLogger;

    /**
     * ParticlePairDpdForceCalculator instance
     */
    private final IParticlePairForceCalculator particlePairDpdForceCalculator;
    
    /**
     * HarmonicBondPropertyCalculator instance
     */
    private final IHarmonicBondForceCalculator harmonicBondConservativeForceCalculator;

    /**
     * ParticlePairElectrostaticsForceConservativeCalculator instance
     */
    private final IParticlePairForceCalculator particlePairElectrostaticsForceConservativeCalculator;
    
    /**
     * True: Force is assigned to f, false: Otherwise
     */
    private final boolean isFAccumulation;
    
    /**
     * True: Force is assigned to ftwo, false: Otherwise
     */
    private final boolean isFtwoAccumulation;
    
    /**
     * Calculation mode for this.particlePairDpdForceCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairDpdForceCalculationMode;

    /**
     * Calculation mode for this.particlePairElectrostaticsForceConservativeCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairElectrostaticsForceConservativeCalculationMode;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation logger
     * @param aParticlePairDpdForceCalculator ParticlePairDpdForceCalculator instance (NOT allowed to be null)
     * @param aParticlePairDpdForceCalculationMode Calculation mode for ParticlePairDpdForceCalculator instance
     * @param aHarmonicBondForceConservativeCalculator HarmonicBondForceConservativeCalculator instance (may be null)
     * @param aParticlePairElectrostaticsForceConservativeCalculator ParticlePairElectrostaticsForceConservativeCalculator instance (may be null)
     * @param aParticlePairElectrostaticsForceConservativeCalculationMode Calculation mode for ParticlePairElectrostaticsForceConservativeCalculator instance
     */
    public ForceAccumulator(
        ILogger aSimulationLogger, 
        IParticlePairForceCalculator aParticlePairDpdForceCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairDpdForceCalculationMode,
        IHarmonicBondForceCalculator aHarmonicBondForceConservativeCalculator,
        IParticlePairForceCalculator aParticlePairElectrostaticsForceConservativeCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairElectrostaticsForceConservativeCalculationMode
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("ForceAccumulator.Constructor: aSimulationLogger is null.");
        }
        if (aParticlePairDpdForceCalculator == null) {
            throw new IllegalArgumentException("ForceAccumulator.Constructor: aParticlePairDpdForceCalculator is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        this.particlePairDpdForceCalculator = aParticlePairDpdForceCalculator;
        this.particlePairDpdForceCalculationMode = aParticlePairDpdForceCalculationMode;
        this.harmonicBondConservativeForceCalculator = aHarmonicBondForceConservativeCalculator;
        this.particlePairElectrostaticsForceConservativeCalculator = aParticlePairElectrostaticsForceConservativeCalculator;
        this.particlePairElectrostaticsForceConservativeCalculationMode = aParticlePairElectrostaticsForceConservativeCalculationMode;
        this.isFAccumulation = 
            this.particlePairDpdForceCalculator.isFAccumulation() || 
            (this.harmonicBondConservativeForceCalculator != null && this.harmonicBondConservativeForceCalculator.isFAccumulation()) ||
            (this.particlePairElectrostaticsForceConservativeCalculator != null && this.particlePairElectrostaticsForceConservativeCalculator.isFAccumulation());
        this.isFtwoAccumulation = 
            this.particlePairDpdForceCalculator.isFtwoAccumulation() || 
            (this.harmonicBondConservativeForceCalculator != null && this.harmonicBondConservativeForceCalculator.isFtwoAccumulation()) ||
            (this.particlePairElectrostaticsForceConservativeCalculator != null && this.particlePairElectrostaticsForceConservativeCalculator.isFtwoAccumulation());
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ForceAccumulator.Constructor");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Accumulates forces to 
     * aParameters.getParticleArrays().getF_x() etc. and
     * aParameters.getParticleArrays().getFtwo_x() etc.
     * (No checks are performed)
     * 
     * @param aBondChunkArraysList List of bond chunk arrays (may be null)
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters
     * @return True: Operation successful, false: Otherwise
     */
    public boolean accumulate_f_and_fTwo(
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters
    ) {
        // <editor-fold defaultstate="collapsed" desc="Initialize f and ftwo if necessary">
        if (this.isFAccumulation) {
            Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
        }
        if (this.isFtwoAccumulation) {
            Arrays.fill(aParameters.getParticleArrays().getFtwo_x(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getFtwo_y(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getFtwo_z(), 0.0);
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="DPD forces">
        boolean tmpIsParticlePairDpdForceCalculationSuccessful = 
            this.particlePairDpdForceCalculator.calculateParticlePairInteractions(
                aR_x,
                aR_y,
                aR_z,
                aParameters,
                this.particlePairDpdForceCalculationMode
            );
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Bond forces">
        boolean tmpIsHarmonicBondForceConservativeCalculationSuccessful = true;
        if (this.harmonicBondConservativeForceCalculator != null) {
            tmpIsHarmonicBondForceConservativeCalculationSuccessful = 
                this.harmonicBondConservativeForceCalculator.calculateBondProperties(
                    aBondChunkArraysList, 
                    aR_x, 
                    aR_y, 
                    aR_z, 
                    aParameters
                );
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Electrostatics forces">
        boolean tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful = true;
        if (this.particlePairElectrostaticsForceConservativeCalculator != null) {
            tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful = 
                this.particlePairElectrostaticsForceConservativeCalculator.calculateParticlePairInteractions(
                    aParameters.getParticleArrays().getChargedParticles_r_x(),
                    aParameters.getParticleArrays().getChargedParticles_r_y(),
                    aParameters.getParticleArrays().getChargedParticles_r_z(),
                    aParameters,
                    this.particlePairElectrostaticsForceConservativeCalculationMode
                );
        }
        // </editor-fold>
        return tmpIsParticlePairDpdForceCalculationSuccessful 
            && tmpIsHarmonicBondForceConservativeCalculationSuccessful
            && tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful;
    }
    
    /**
     * Executor service shutdown
     */
    public void shutdownExecutorService() {
        this.particlePairDpdForceCalculator.shutdownExecutorService();
        if (this.harmonicBondConservativeForceCalculator != null) {
            this.harmonicBondConservativeForceCalculator.shutdownExecutorService();
        }
        if (this.particlePairElectrostaticsForceConservativeCalculator != null) {
            this.particlePairElectrostaticsForceConservativeCalculator.shutdownExecutorService();
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * ParticlePairDpdForceCalculator instance
     * 
     * @return ParticlePairDpdForceCalculator instance
     */
    public IParticlePairForceCalculator getParticlePairDpdForceCalculator() {
        return this.particlePairDpdForceCalculator;
    }
    
    /**
     * HarmonicBondForceConservativeCalculator instance
     * 
     * @return HarmonicBondForceConservativeCalculator instance
     */
    public IHarmonicBondForceCalculator getHarmonicBondForceConservativeCalculator() {
        return this.harmonicBondConservativeForceCalculator;
    }

    /**
     * ParticlePairElectrostaticsForceConservativeCalculator instance
     * 
     * @return ParticlePairElectrostaticsForceConservativeCalculator instance
     */
    public IParticlePairForceCalculator getParticlePairElectrostaticsForceConservativeCalculator() {
        return this.particlePairElectrostaticsForceConservativeCalculator;
    }
    // </editor-fold>
    
}
