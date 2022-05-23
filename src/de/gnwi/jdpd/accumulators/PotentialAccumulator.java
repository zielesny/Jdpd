/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2022  Achim Zielesny (achim.zielesny@googlemail.com)
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
import de.gnwi.jdpd.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.ExtendedAdderGroup;
import java.util.LinkedList;

/**
 * Potential accumulator that uses ExtendedAdderGroup instance for 
 * accumulative ("+=") operations.
 * 
 * @author Achim Zielesny
 */
public class PotentialAccumulator {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Simulation logger
     */
    private final ILogger simulationLogger;

    /**
     * ParticlePairDpdPotentialCalculator instance
     */
    private final IParticlePairInteractionCalculator particlePairDpdPotentialCalculator;
    
    /**
     * HarmonicBondPotentialCalculator instance
     */
    private final IHarmonicBondPropertyCalculator harmonicBondPotentialCalculator;

    /**
     * ParticlePairElectrostaticsPotentialCalculator instance
     */
    private final IParticlePairInteractionCalculator particlePairElectrostaticsPotentialCalculator;
    
    /**
     * Extended adder group
     */
    private final ExtendedAdderGroup extendedAdderGroup;
    
    /**
     * Calculation mode for this.particlePairDpdForceCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairDpdPotentialCalculationMode;

    /**
     * Calculation mode for this.particlePairElectrostaticsForceConservativeCalculator
     */
    private final ParticlePairInteractionCalculator.CellBasedCalculationMode particlePairElectrostaticsPotentialCalculationMode;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation logger
     * @param aParticlePairDpdPotentialCalculator ParticlePairDpdPotentialCalculator instance (NOT allowed to be null)
     * @param aParticlePairDpdPotentialCalculationMode Calculation mode for ParticlePairDpdPotentialCalculator instance
     * @param aHarmonicBondPotentialCalculator HarmonicBondPotentialCalculator instance (may be null)
     * @param aParticlePairElectrostaticsPotentialCalculator ParticlePairElectrostaticsPotentialCalculator instance (may be null)
     * @param aParticlePairElectrostaticsPotentialCalculationMode Calculation mode for aParticlePairElectrostaticsPotentialCalculator instance
     */
    public PotentialAccumulator(
        ILogger aSimulationLogger, 
        IParticlePairInteractionCalculator aParticlePairDpdPotentialCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairDpdPotentialCalculationMode,
        IHarmonicBondPropertyCalculator aHarmonicBondPotentialCalculator,
        IParticlePairInteractionCalculator aParticlePairElectrostaticsPotentialCalculator,
        ParticlePairInteractionCalculator.CellBasedCalculationMode aParticlePairElectrostaticsPotentialCalculationMode
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("PotentialAccumulator.Constructor: aSimulationLogger is null.");
        }
        if (aParticlePairDpdPotentialCalculator == null) {
            throw new IllegalArgumentException("PotentialAccumulator.Constructor: aParticlePairDpdPotentialCalculator is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        this.particlePairDpdPotentialCalculator = aParticlePairDpdPotentialCalculator;
        this.particlePairDpdPotentialCalculationMode = aParticlePairDpdPotentialCalculationMode;
        this.harmonicBondPotentialCalculator = aHarmonicBondPotentialCalculator;
        this.particlePairElectrostaticsPotentialCalculator = aParticlePairElectrostaticsPotentialCalculator;
        this.particlePairElectrostaticsPotentialCalculationMode = aParticlePairElectrostaticsPotentialCalculationMode;
        this.extendedAdderGroup = new ExtendedAdderGroup();
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("PotentialAccumulator.Constructor");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Accumulates potentials to 
     * aParameters.getParticleArrays().getCalculationAdders()
     * (No checks are performed)
     * 
     * @param aBondChunkArraysList List of bond chunk arrays (may be null)
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters
     * @return True: Operation successful, false: Otherwise
     */
    public boolean accumulatePotentials(
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters
    ) {
        // <editor-fold defaultstate="collapsed" desc="Reset extend adder group">
        this.extendedAdderGroup.reset();
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="DPD potentials">
        boolean tmpIsParticlePairDpdPotentialCalculationSuccessful = 
            this.particlePairDpdPotentialCalculator.calculateParticlePairInteractions(
                aR_x,
                aR_y,
                aR_z,
                aParameters,
                this.particlePairDpdPotentialCalculationMode
            );
        double tmpUpotDpd = this.particlePairDpdPotentialCalculator.getAccumulatedPotentialEnergyAddersSum();
        double tmpPressureXDpd = this.particlePairDpdPotentialCalculator.getAccumulatedPressureXAddersSum();
        double tmpPressureYDpd = this.particlePairDpdPotentialCalculator.getAccumulatedPressureYAddersSum();
        double tmpPressureZDpd = this.particlePairDpdPotentialCalculator.getAccumulatedPressureZAddersSum();
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: uPot(DPD)      = " + String.valueOf(tmpUpotDpd));
        this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureX(DPD) = " + String.valueOf(tmpPressureXDpd));
        this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureY(DPD) = " + String.valueOf(tmpPressureYDpd));
        this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureZ(DPD) = " + String.valueOf(tmpPressureZDpd));
        // </editor-fold>
        this.extendedAdderGroup.getPotentialEnergyAdder().add(tmpUpotDpd);
        this.extendedAdderGroup.getDpdPotentialEnergyAdder().add(tmpUpotDpd);
        this.extendedAdderGroup.getPressureXAdder().add(tmpPressureXDpd);
        this.extendedAdderGroup.getPressureYAdder().add(tmpPressureYDpd);
        this.extendedAdderGroup.getPressureZAdder().add(tmpPressureZDpd);
        this.extendedAdderGroup.getDpdPressureXAdder().add(tmpPressureXDpd);
        this.extendedAdderGroup.getDpdPressureYAdder().add(tmpPressureYDpd);
        this.extendedAdderGroup.getDpdPressureZAdder().add(tmpPressureZDpd);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Bond potentials">
        boolean tmpIsHarmonicBondPotentialCalculationSuccessful = true;
        if (this.harmonicBondPotentialCalculator != null) {
            tmpIsHarmonicBondPotentialCalculationSuccessful = 
                this.harmonicBondPotentialCalculator.calculateBondProperties(
                    aBondChunkArraysList, 
                    aR_x, 
                    aR_y, 
                    aR_z, 
                    aParameters
                );
            double tmpUpotBonds = this.harmonicBondPotentialCalculator.getAccumulatedPotentialEnergyAddersSum();
            double tmpPressureXBonds = this.harmonicBondPotentialCalculator.getAccumulatedPressureXAddersSum();
            double tmpPressureYBonds = this.harmonicBondPotentialCalculator.getAccumulatedPressureYAddersSum();
            double tmpPressureZBonds = this.harmonicBondPotentialCalculator.getAccumulatedPressureZAddersSum();
            // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: uPot(Bonds)      = " + String.valueOf(tmpUpotBonds));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureX(Bonds) = " + String.valueOf(tmpPressureXBonds));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureY(Bonds) = " + String.valueOf(tmpPressureYBonds));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureZ(Bonds) = " + String.valueOf(tmpPressureZBonds));
            // </editor-fold>
            this.extendedAdderGroup.getPotentialEnergyAdder().add(tmpUpotBonds);
            this.extendedAdderGroup.getBondPotentialEnergyAdder().add(tmpUpotBonds);
            this.extendedAdderGroup.getPressureXAdder().add(tmpPressureXBonds);
            this.extendedAdderGroup.getPressureYAdder().add(tmpPressureYBonds);
            this.extendedAdderGroup.getPressureZAdder().add(tmpPressureZBonds);
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Electrostatics potentials">
        boolean tmpIsParticlePairElectrostaticsPotentialCalculationSuccessful = true;
        if (this.particlePairElectrostaticsPotentialCalculator != null) {
            tmpIsParticlePairElectrostaticsPotentialCalculationSuccessful = 
                this.particlePairElectrostaticsPotentialCalculator.calculateParticlePairInteractions(
                    aParameters.getParticleArrays().getChargedParticles_r_x(),
                    aParameters.getParticleArrays().getChargedParticles_r_y(),
                    aParameters.getParticleArrays().getChargedParticles_r_z(),
                    aParameters,
                    this.particlePairElectrostaticsPotentialCalculationMode
                );
            double tmpUpotParticlePairElectrostatics = this.particlePairElectrostaticsPotentialCalculator.getAccumulatedPotentialEnergyAddersSum();
            double tmpPressureXParticlePairElectrostatics = this.particlePairElectrostaticsPotentialCalculator.getAccumulatedPressureXAddersSum();
            double tmpPressureYParticlePairElectrostatics = this.particlePairElectrostaticsPotentialCalculator.getAccumulatedPressureYAddersSum();
            double tmpPressureZParticlePairElectrostatics = this.particlePairElectrostaticsPotentialCalculator.getAccumulatedPressureZAddersSum();
            // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: uPot(Estat)      = " + String.valueOf(tmpUpotParticlePairElectrostatics));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureX(Estat) = " + String.valueOf(tmpPressureXParticlePairElectrostatics));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureY(Estat) = " + String.valueOf(tmpPressureYParticlePairElectrostatics));
            this.simulationLogger.appendIntermediateResults("PotentialAccumulator.accumulate: PressureZ(Estat) = " + String.valueOf(tmpPressureZParticlePairElectrostatics));
            // </editor-fold>
            this.extendedAdderGroup.getPotentialEnergyAdder().add(tmpUpotParticlePairElectrostatics);
            this.extendedAdderGroup.getElectrostaticsPotentialEnergyAdder().add(tmpUpotParticlePairElectrostatics);
            this.extendedAdderGroup.getPressureXAdder().add(tmpPressureXParticlePairElectrostatics);
            this.extendedAdderGroup.getPressureYAdder().add(tmpPressureYParticlePairElectrostatics);
            this.extendedAdderGroup.getPressureZAdder().add(tmpPressureZParticlePairElectrostatics);
        }
        // </editor-fold>
        return tmpIsParticlePairDpdPotentialCalculationSuccessful 
            && tmpIsHarmonicBondPotentialCalculationSuccessful
            && tmpIsParticlePairElectrostaticsPotentialCalculationSuccessful;
    }

    /**
     * Executor service shutdown
     */
    public void shutdownExecutorService() {
        this.particlePairDpdPotentialCalculator.shutdownExecutorService();
        if (this.harmonicBondPotentialCalculator != null) {
            this.harmonicBondPotentialCalculator.shutdownExecutorService();
        }
        if (this.particlePairElectrostaticsPotentialCalculator != null) {
            this.particlePairElectrostaticsPotentialCalculator.shutdownExecutorService();
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * ParticlePairDpdPotentialCalculator instance
     * 
     * @return ParticlePairDpdPotentialCalculator instance
     */
    public IParticlePairInteractionCalculator getParticlePairDpdPotentialCalculator() {
        return this.particlePairDpdPotentialCalculator;
    }
    
    /**
     * HarmonicBondPotentialCalculator instance
     * 
     * @return HarmonicBondPotentialCalculator instance
     */
    public IHarmonicBondPropertyCalculator getHarmonicBondPotentialCalculator() {
        return this.harmonicBondPotentialCalculator;
    }

    /**
     * ParticlePairElectrostaticsPotentialCalculator instance
     * 
     * @return ParticlePairElectrostaticsPotentialCalculator instance
     */
    public IParticlePairInteractionCalculator getParticlePairElectrostaticsPotentialCalculator() {
        return this.particlePairElectrostaticsPotentialCalculator;
    }
    
    /**
     * Extended adder group
     * 
     * @return Extended adder group
     */
    public ExtendedAdderGroup getExtendedAdderGroup() {
        return this.extendedAdderGroup;
    }
    // </editor-fold>
    
}
