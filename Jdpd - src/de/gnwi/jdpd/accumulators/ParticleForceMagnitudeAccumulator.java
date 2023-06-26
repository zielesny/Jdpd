/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpd.accumulators;

import de.gnwi.jdpd.samples.harmonicBonds.HarmonicBondChunkArrays;
import de.gnwi.jdpd.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.Utils;
import java.util.Arrays;
import java.util.LinkedList;

/**
 * Accumulator for particle force magnitude analysis.
 * The particle related arrays
 * ParticleArrays.getF_x
 * ParticleArrays.getF_y
 * ParticleArrays.getF_z
 * are used for accumulative ("+=" and "-=") operations but their initial
 * value is saved and finally restored, i.e. they are NOT changed by method 
 * analyzeParticleForceMagnitudes().
 * 
 * @author Achim Zielesny
 */
public class ParticleForceMagnitudeAccumulator {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Simulation logger
     */
    private final ILogger simulationLogger;

    /**
     * ParticlePairDpdForceConservativeCalculator instance
     */
    private final IParticlePairForceCalculator particlePairDpdForceConservativeCalculator;

    /**
     * ParticlePairDpdForceRandomCalculator instance
     */
    private final IParticlePairForceCalculator particlePairDpdForceRandomCalculator;

    /**
     * ParticlePairDpdForceDissipativeCalculator instance
     */
    private final IParticlePairForceCalculator particlePairDpdForceDissipativeCalculator;
    
    /**
     * HarmonicBondPropertyCalculator instance
     */
    private final IHarmonicBondForceCalculator harmonicBondConservativeForceCalculator;

    /**
     * ParticlePairElectrostaticsForceConservativeCalculator instance
     */
    private final IParticlePairForceCalculator particlePairElectrostaticsForceConservativeCalculator;
    
    /**
     * Calculation mode for this.particlePairDpdForceConservativeCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairDpdForceConservativeCalculationMode;

    /**
     * Calculation mode for this.particlePairElectrostaticsForceConservativeCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairElectrostaticsForceConservativeCalculationMode;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Auxiliary force x-component
     */
    private double[] fAux_x;

    /**
     * Auxiliary force y-component
     */
    private double[] fAux_y;
    
    /**
     * Auxiliary force z-component
     */
    private double[] fAux_z;

    /**
     * Auxiliary total force x-component
     */
    private double[] fAuxTot_x;

    /**
     * Auxiliary total force y-component
     */
    private double[] fAuxTot_y;
    
    /**
     * Auxiliary total force z-component
     */
    private double[] fAuxTot_z;
    
    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative DPD force magnitudes of particles
     */
    private double[] minMeanMaxDpdForceConservativeMagnitude;

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of random DPD force magnitudes of particles
     */
    private double[] minMeanMaxDpdForceRandomMagnitude;

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of dissipative DPD force magnitudes of particles
     */
    private double[] minMeanMaxDpdForceDissipativeMagnitude;

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative bond force magnitudes of 
     * particles
     */
    private double[] minMeanMaxBondForceMagnitude;

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of electrostatics conservative force magnitudes 
     * of particles
     */
    private double[] minMeanMaxElectrostaticsForceMagnitude;

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of total force magnitudes 
     * of particles
     */
    private double[] minMeanMaxTotalForceMagnitude;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation logger
     * @param aParticlePairDpdForceConservativeCalculator ParticlePairDpdForceConservativeCalculator instance (NOT allowed to be null)
     * @param aParticlePairDpdForceConservativeCalculationMode Calculation mode for ParticlePairDpdForceConservativeCalculator instance
     * @param aParticlePairDpdForceRandomCalculator ParticlePairDpdForceRandomCalculator instance (NOT allowed to be null)
     * @param aParticlePairDpdForceDissipativeCalculator ParticlePairDpdForceDissipativeCalculator instance (NOT allowed to be null)
     * @param aHarmonicBondForceConservativeCalculator HarmonicBondForceConservativeCalculator instance (may be null)
     * @param aParticlePairElectrostaticsForceConservativeCalculator ParticlePairElectrostaticsForceConservativeCalculator instance (may be null)
     * @param aParticlePairElectrostaticsForceConservativeCalculationMode Calculation mode for ParticlePairElectrostaticsForceConservativeCalculator instance
     */
    public ParticleForceMagnitudeAccumulator(
        ILogger aSimulationLogger, 
        IParticlePairForceCalculator aParticlePairDpdForceConservativeCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairDpdForceConservativeCalculationMode,
        IParticlePairForceCalculator aParticlePairDpdForceRandomCalculator,
        IParticlePairForceCalculator aParticlePairDpdForceDissipativeCalculator,
        IHarmonicBondForceCalculator aHarmonicBondForceConservativeCalculator,
        IParticlePairForceCalculator aParticlePairElectrostaticsForceConservativeCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairElectrostaticsForceConservativeCalculationMode
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: aSimulationLogger is null.");
        }
        if (aParticlePairDpdForceConservativeCalculator == null) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: aParticlePairDpdForceConservativeCalculator is null.");
        }
        if (!aParticlePairDpdForceConservativeCalculator.getParticlePairDistanceParametersCacheActivity()) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: Cache of aParticlePairDpdForceConservativeCalculator is not active.");
        }
        if (aParticlePairDpdForceRandomCalculator == null) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: aParticlePairDpdForceRandomCalculator is null.");
        }
        if (aParticlePairDpdForceDissipativeCalculator == null) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: aParticlePairDpdForceDissipativeCalculator is null.");
        }
        boolean tmpIsFAccumulation = 
            aParticlePairDpdForceConservativeCalculator.isFAccumulation() || 
            (aHarmonicBondForceConservativeCalculator != null && aHarmonicBondForceConservativeCalculator.isFAccumulation()) ||
            (aParticlePairElectrostaticsForceConservativeCalculator != null && aParticlePairElectrostaticsForceConservativeCalculator.isFAccumulation());
        if (!tmpIsFAccumulation) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: F accumulation MUST be supported.");
        }
        boolean tmpIsFtwoAccumulation = 
            aParticlePairDpdForceConservativeCalculator.isFtwoAccumulation() || 
            (aHarmonicBondForceConservativeCalculator != null && aHarmonicBondForceConservativeCalculator.isFtwoAccumulation()) ||
            (aParticlePairElectrostaticsForceConservativeCalculator != null && aParticlePairElectrostaticsForceConservativeCalculator.isFtwoAccumulation());
        if (tmpIsFtwoAccumulation) {
            throw new IllegalArgumentException("ParticleForceMagnitudeAccumulator.Constructor: Ftwo accumulation is NOT supported.");
        }
        // </editor-fold>
        
        this.simulationLogger = aSimulationLogger;
        this.particlePairDpdForceConservativeCalculator = aParticlePairDpdForceConservativeCalculator;
        this.particlePairDpdForceConservativeCalculationMode = aParticlePairDpdForceConservativeCalculationMode;
        this.particlePairDpdForceRandomCalculator = aParticlePairDpdForceRandomCalculator;
        this.particlePairDpdForceDissipativeCalculator = aParticlePairDpdForceDissipativeCalculator;
        this.harmonicBondConservativeForceCalculator = aHarmonicBondForceConservativeCalculator;
        this.particlePairElectrostaticsForceConservativeCalculator = aParticlePairElectrostaticsForceConservativeCalculator;
        this.particlePairElectrostaticsForceConservativeCalculationMode = aParticlePairElectrostaticsForceConservativeCalculationMode;

        this.minMeanMaxDpdForceConservativeMagnitude = null;
        this.minMeanMaxDpdForceRandomMagnitude = null;
        this.minMeanMaxDpdForceDissipativeMagnitude = null;
        this.minMeanMaxBondForceMagnitude = null;
        this.minMeanMaxElectrostaticsForceMagnitude = null;
        this.minMeanMaxTotalForceMagnitude = null;

        this.fAux_x = null;
        this.fAux_y = null;
        this.fAux_z = null;
        this.fAuxTot_x = null;
        this.fAuxTot_y = null;
        this.fAuxTot_z = null;
        
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticleForceMagnitudeAccumulator.Constructor");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Analyses particle force magnitudes
     * NOTE: The particle related arrays
     * ParticleArrays.getF_x
     * ParticleArrays.getF_y
     * ParticleArrays.getF_z
     * are used for accumulative ("+=" and "-=") operations but their initial
     * value is saved and finally restored, i.e. they are NOT changed.
     * (No checks are performed)
     * 
     * @param aBondChunkArraysList List of bond chunk arrays (may be null)
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters
     * @return True: Operation successful, false: Otherwise
     */
    public boolean analyzeParticleForceMagnitudes(
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters
    ) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Save ParticleArrays.getF_x/y/z">
            // Save ParticleArrays.getF_x/y/z in this.fAux_x/y/z
            if (this.fAux_x == null) {
                this.fAux_x = new double[aParameters.getParticleArrays().getF_x().length];
                this.fAux_y = new double[aParameters.getParticleArrays().getF_x().length];
                this.fAux_z = new double[aParameters.getParticleArrays().getF_x().length];
                this.fAuxTot_x = new double[aParameters.getParticleArrays().getF_x().length];
                this.fAuxTot_y = new double[aParameters.getParticleArrays().getF_x().length];
                this.fAuxTot_z = new double[aParameters.getParticleArrays().getF_x().length];
            }
            Arrays.fill(this.fAuxTot_x, 0.0);
            Arrays.fill(this.fAuxTot_y, 0.0);
            Arrays.fill(this.fAuxTot_z, 0.0);
            Utils.copyToOld(
                aParameters.getParticleArrays().getF_x(), 
                aParameters.getParticleArrays().getF_y(), 
                aParameters.getParticleArrays().getF_z(), 
                this.fAux_x, 
                this.fAux_y, 
                this.fAux_z
            );
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Conservative DPD forces">
            // Initialize ParticleArrays.getF_x/y/z
            Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
            boolean tmpIsParticlePairDpdForceConservativeCalculationSuccessful = 
                this.particlePairDpdForceConservativeCalculator.calculateParticlePairInteractions(
                    aR_x,
                    aR_y,
                    aR_z,
                    aParameters,
                    this.particlePairDpdForceConservativeCalculationMode
                );
            this.minMeanMaxDpdForceConservativeMagnitude = 
                Utils.calculateMinMeanMax(
                    aParameters.getParticleArrays().getF_x(),
                    aParameters.getParticleArrays().getF_y(),
                    aParameters.getParticleArrays().getF_z()
                );
            Utils.add_b_to_a(
                this.fAuxTot_x,
                this.fAuxTot_y,
                this.fAuxTot_z,
                aParameters.getParticleArrays().getF_x(),
                aParameters.getParticleArrays().getF_y(),
                aParameters.getParticleArrays().getF_z()
            );
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Random DPD forces (with cache)">
            // Initialize ParticleArrays.getF_x/y/z
            Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
            this.particlePairDpdForceRandomCalculator.setParticlePairDistanceParametersCache(this.particlePairDpdForceConservativeCalculator);
            boolean tmpIsParticlePairDpdForceRandomCalculationSuccessful = 
                this.particlePairDpdForceRandomCalculator.calculateParticlePairInteractions(
                    aR_x,
                    aR_y,
                    aR_z,
                    aParameters,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_CACHE
                );
            this.minMeanMaxDpdForceRandomMagnitude = 
                Utils.calculateMinMeanMax(
                    aParameters.getParticleArrays().getF_x(),
                    aParameters.getParticleArrays().getF_y(),
                    aParameters.getParticleArrays().getF_z()
                );
            Utils.add_b_to_a(
                this.fAuxTot_x,
                this.fAuxTot_y,
                this.fAuxTot_z,
                aParameters.getParticleArrays().getF_x(),
                aParameters.getParticleArrays().getF_y(),
                aParameters.getParticleArrays().getF_z()
            );
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Dissipative DPD forces (with cache)">
            // Initialize ParticleArrays.getF_x/y/z
            Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
            Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
            this.particlePairDpdForceDissipativeCalculator.setParticlePairDistanceParametersCache(this.particlePairDpdForceConservativeCalculator);
            boolean tmpIsParticlePairDpdForceDissipativeCalculationSuccessful = 
                this.particlePairDpdForceDissipativeCalculator.calculateParticlePairInteractions(
                    aR_x,
                    aR_y,
                    aR_z,
                    aParameters,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_CACHE
                );
            this.minMeanMaxDpdForceDissipativeMagnitude = 
                Utils.calculateMinMeanMax(
                    aParameters.getParticleArrays().getF_x(),
                    aParameters.getParticleArrays().getF_y(),
                    aParameters.getParticleArrays().getF_z()
                );
            Utils.add_b_to_a(
                this.fAuxTot_x,
                this.fAuxTot_y,
                this.fAuxTot_z,
                aParameters.getParticleArrays().getF_x(),
                aParameters.getParticleArrays().getF_y(),
                aParameters.getParticleArrays().getF_z()
            );
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Bond forces">
            boolean tmpIsHarmonicBondForceConservativeCalculationSuccessful = true;
            if (this.harmonicBondConservativeForceCalculator != null) {
                Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
                Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
                Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
                tmpIsHarmonicBondForceConservativeCalculationSuccessful = 
                    this.harmonicBondConservativeForceCalculator.calculateBondProperties(
                        aBondChunkArraysList, 
                        aR_x, 
                        aR_y, 
                        aR_z, 
                        aParameters
                    );
                this.minMeanMaxBondForceMagnitude = 
                    Utils.calculateMinMeanMax(
                        aParameters.getParticleArrays().getF_x(),
                        aParameters.getParticleArrays().getF_y(),
                        aParameters.getParticleArrays().getF_z()
                    );
                Utils.add_b_to_a(
                    this.fAuxTot_x,
                    this.fAuxTot_y,
                    this.fAuxTot_z,
                    aParameters.getParticleArrays().getF_x(),
                    aParameters.getParticleArrays().getF_y(),
                    aParameters.getParticleArrays().getF_z()
                );
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Electrostatics forces">
            boolean tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful = true;
            if (this.particlePairElectrostaticsForceConservativeCalculator != null) {
                Arrays.fill(aParameters.getParticleArrays().getF_x(), 0.0);
                Arrays.fill(aParameters.getParticleArrays().getF_y(), 0.0);
                Arrays.fill(aParameters.getParticleArrays().getF_z(), 0.0);
                tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful = 
                    this.particlePairElectrostaticsForceConservativeCalculator.calculateParticlePairInteractions(
                        aParameters.getParticleArrays().getChargedParticles_r_x(),
                        aParameters.getParticleArrays().getChargedParticles_r_y(),
                        aParameters.getParticleArrays().getChargedParticles_r_z(),
                        aParameters,
                        this.particlePairElectrostaticsForceConservativeCalculationMode
                    );
                this.minMeanMaxElectrostaticsForceMagnitude = 
                    Utils.calculateMinMeanMax(
                        aParameters.getParticleArrays().getF_x(),
                        aParameters.getParticleArrays().getF_y(),
                        aParameters.getParticleArrays().getF_z()
                    );
                Utils.add_b_to_a(
                    this.fAuxTot_x,
                    this.fAuxTot_y,
                    this.fAuxTot_z,
                    aParameters.getParticleArrays().getF_x(),
                    aParameters.getParticleArrays().getF_y(),
                    aParameters.getParticleArrays().getF_z()
                );
            }
            this.minMeanMaxTotalForceMagnitude = 
                Utils.calculateMinMeanMax(
                    this.fAuxTot_x,
                    this.fAuxTot_y,
                    this.fAuxTot_z
                );
            // </editor-fold>
            return tmpIsParticlePairDpdForceConservativeCalculationSuccessful 
                && tmpIsParticlePairDpdForceRandomCalculationSuccessful
                && tmpIsParticlePairDpdForceDissipativeCalculationSuccessful
                && tmpIsHarmonicBondForceConservativeCalculationSuccessful
                && tmpIsParticlePairElectrostaticsForceConservativeCalculationSuccessful;
        } finally {
            // <editor-fold defaultstate="collapsed" desc="Restore ParticleArrays.getF_x/y/z">
            // Restore ParticleArrays.getF_x/y/z from this.fAux_x/y/z
            Utils.copyToOld(
                this.fAux_x, 
                this.fAux_y, 
                this.fAux_z,
                aParameters.getParticleArrays().getF_x(), 
                aParameters.getParticleArrays().getF_y(), 
                aParameters.getParticleArrays().getF_z()
            );
            // </editor-fold>
        }
    }
    
    /**
     * Executor service shutdown
     */
    public void shutdownExecutorService() {
        this.particlePairDpdForceConservativeCalculator.shutdownExecutorService();
        this.particlePairDpdForceRandomCalculator.shutdownExecutorService();
        this.particlePairDpdForceDissipativeCalculator.shutdownExecutorService();
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
     * ParticlePairDpdForceConservativeCalculator instance
     * 
     * @return ParticlePairDpdForceConservativeCalculator instance
     */
    public IParticlePairForceCalculator getParticlePairDpdForceConservativeCalculator() {
        return this.particlePairDpdForceConservativeCalculator;
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

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative DPD force magnitudes of particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative DPD force magnitudes of particles
     */
    public double[] getMinMeanMaxDpdForceConservativeMagnitude() {
        return this.minMeanMaxDpdForceConservativeMagnitude;
    }

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of random DPD force magnitudes of particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of random DPD force magnitudes of particles
     */
    public double[] getMinMeanMaxDpdForceRandomMagnitude() {
        return this.minMeanMaxDpdForceRandomMagnitude;
    }

    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of dissipative DPD force magnitudes of particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of dissipative DPD force magnitudes of particles
     */
    public double[] getMinMeanMaxDpdForceDissipativeMagnitude() {
        return this.minMeanMaxDpdForceDissipativeMagnitude;
    }
    
    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative bond force magnitudes of 
     * particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of conservative bond force magnitudes of 
     * particles
     */
    public double[] getMinMeanMaxBondForceMagnitude() {
        return this.minMeanMaxBondForceMagnitude;
    }
    
    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of electrostatics conservative force magnitudes 
     * of particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of electrostatics conservative force magnitudes 
     * of particles
     */
    public double[] getMinMeanMaxElectrostaticsForceMagnitude() {
        return this.minMeanMaxElectrostaticsForceMagnitude;
    }
    
    /**
     * Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of total force magnitudes 
     * of particles
     * 
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of total force magnitudes 
     * of particles
     */
    public double[] getMinMeanMaxTotalForceMagnitude() {
        return this.minMeanMaxTotalForceMagnitude;
    }
    // </editor-fold>
    
}
