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
package de.gnwi.jdpdsp;

import de.gnwi.jdpd.interfaces.IDpdSimulationTask;
import de.gnwi.jdpd.interfaces.ILogLevel;
import de.gnwi.jdpd.interfaces.IProgressMonitor;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.utilities.Factory;
import de.gnwi.jdpdsp.samples.harmonicBonds.ParticlePairHarmonicBond;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBond;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondChunkArrays;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondChunkGenerator;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.utilities.Utils;
import de.gnwi.jdpdsp.interfaces.IInput;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IOutput;
import de.gnwi.jdpdsp.parameters.ParticleTypes;
import de.gnwi.jdpdsp.parameters.ParticleArrays;
import de.gnwi.jdpdsp.parameters.ChemicalSystemDescription;
import de.gnwi.jdpdsp.parameters.InteractionDescription;
import de.gnwi.jdpdsp.parameters.MoleculeTypes;
import de.gnwi.jdpdsp.parameters.RestartInfo;
import de.gnwi.jdpdsp.parameters.SimulationCounts;
import de.gnwi.jdpdsp.parameters.SimulationDescription;
import de.gnwi.jdpdsp.parameters.TestObjects;
import de.gnwi.jdpdsp.utilities.MoleculeDescription;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.atomic.AtomicInteger;
import de.gnwi.jdpdsp.interfaces.IRandom;
import de.gnwi.jdpdsp.particlePosition.ParticlePositionPool;
import java.util.concurrent.Callable;
import de.gnwi.jdpdsp.interfaces.ITimeStepCalculator;
import de.gnwi.jdpdsp.movement.MoleculeFixationInfo;
import de.gnwi.jdpdsp.movement.MoleculeVelocityFixationInfo;
import de.gnwi.jdpdsp.rg.MoleculeRgValue;
import de.gnwi.jdpdsp.rg.RgCalculator;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.movement.MoleculeFixationDescription;
import de.gnwi.jdpdsp.rg.MoleculeRgCalculationDescription;
import de.gnwi.jdpdsp.movement.MoleculeVelocityFixationDescription;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.movement.MoleculeAccelerationDescription;
import de.gnwi.jdpdsp.movement.MoleculeAccelerationInfo;
import de.gnwi.jdpdsp.movement.MoleculeBoundaryDescription;
import de.gnwi.jdpdsp.movement.MoleculeBoundaryInfo;
import de.gnwi.jdpdsp.movement.MoleculeSphereDescription;
import de.gnwi.jdpdsp.movement.MoleculeSphereInfo;
import de.gnwi.jdpdsp.nearestNeighbor.NearestNeighborManager;
import de.gnwi.jdpdsp.nearestNeighbor.NearestNeighborBaseParticleDescription;
import de.gnwi.jdpdsp.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.utilities.Electrostatics;
import java.util.Arrays;
import org.apache.commons.math3.util.FastMath;

/**
 * Task for Molecular Fragment DPD simulation performed with single 
 * precision (SP) arithmetic.
 *
 * @author Achim Zielesny
 *
 */
public class DpdSimulationTaskSP implements Callable<Boolean>, IDpdSimulationTask {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Numeric constants
     */
    private static final float ONE_A_HALF = 1.5f;
    private static final float HUNDRED = 100.0f;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Parallelisation info
     */
    private final ParallelizationInfo parallelizationInfo;

    /**
     * Maximum number of position correction trials
     */
    private final int maximumNumberOfPositionCorrectionTrials;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Factory for new objects
     */
    private Factory factory;
    
    /**
     * Simulation output
     */
    private IOutput simulationOutput;
    
    /**
     * Simulation logger
     */
    private ILogger simulationLogger;
    
    /**
     * Simulation progress monitor
     */
    private IProgressMonitor progressMonitor;
    
    /** 
     * Simulation parameters
     */
    private Parameters parameters;
    
    /**
     * Random number seed
     */
    private AtomicInteger randomNumberSeed;
    
    /**
     * True: Simulation is stopped, false: Otherwise
     */
    private boolean isSimulationStopped;
    
    /**
     * Particle position pools
     */
    private ParticlePositionPool particlePositionPool;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     *
     * @param aRestartInfo Restart info (may be null)
     * @param aSimulationInput Simulation input
     * @param aSimulationOutput Simulation output
     * @param aProgressMonitor Simulation progress monitor
     * @param aSimulationLogger Simulation logger
     * @param aParallelizationInfo Parallelisation info
     * @throws IllegalArgumentException Thrown if argument is invalid
     * @throws Exception Thrown if exception occurred
     */
    public DpdSimulationTaskSP(
        RestartInfo aRestartInfo,
        IInput aSimulationInput, 
        IOutput aSimulationOutput,
        IProgressMonitor aProgressMonitor,
        ILogger aSimulationLogger,
        ParallelizationInfo aParallelizationInfo
    ) throws IllegalArgumentException, Exception {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationInput == null) {
            throw new IllegalArgumentException("DpdSimulationTask.Constructor: aSimulationInput is null.");
        }
        if (aSimulationOutput == null) {
            throw new IllegalArgumentException("DpdSimulationTask.Constructor: aSimulationOutput is null.");
        }
        if (aProgressMonitor == null) {
            throw new IllegalArgumentException("DpdSimulationTask.Constructor: aProgressMonitor is null.");
        }
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("DpdSimulationTask.Constructor: aSimulationLogger is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("DpdSimulationTask.Constructor: aParallelizationInfo is null.");
        }
        // </editor-fold>
        // FIRST: Set and start simulation logger
        this.simulationLogger = aSimulationLogger;
        this.simulationLogger.start();
        this.parallelizationInfo = aParallelizationInfo;
        this.factory = aSimulationInput.getFactory();
        this.simulationOutput = aSimulationOutput;
        this.progressMonitor = aProgressMonitor;
        try {
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            long tmpId = this.simulationLogger.getId();
            this.simulationLogger.appendMethodCallStart("DpdSimulationTask.Constructor", tmpId);
            // </editor-fold>
            // 1. Initialise class variables
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.NOT_STARTED);
            this.isSimulationStopped = false;
            this.particlePositionPool = new ParticlePositionPool();
            this.maximumNumberOfPositionCorrectionTrials = aSimulationInput.getMaximumNumberOfPositionCorrectionTrials();
            // IMPORTANT: Set particle position pool in simulation output for re-use operation after output
            this.simulationOutput.setParticlePositionPool(this.particlePositionPool);
            // 2. Initialise simulation parameters
            this.setParameters(aSimulationInput, aRestartInfo);
            // 3. Correct molecule boundaries if necessary
            this.correctMoleculeBoundaries();
            // 4. Correct molecule spheres if necessary
            this.correctMoleculeSpheres();
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            this.simulationLogger.appendMethodCallEnd("DpdSimulationTask.Constructor", tmpId);
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("DpdSimulationTask.Constructor", Utils.getStacktrace(anException));
            // </editor-fold>
            this.releaseMemoryAndFinishLogging();
            throw anException;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Stops simulation
     */
    @Override
    public void stopSimulation() {
        this.isSimulationStopped = true;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Overriden call() methods">
    /**
     * Executes specified Molecular-Fragment DPD simulation
     *
     * @return True: Simulation has been successfully executed, false: Otherwise.
     * @throws Exception Thrown if error occurred
     */
    @Override
    public Boolean call() throws Exception {
        try {
            // <editor-fold defaultstate="collapsed" desc="Simulator starts: Set simulator progress to 0 percent">
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.STARTED);
            this.progressMonitor.setProgressInPercent(0);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            long tmpIdDoInBackground = this.simulationLogger.getId();
            this.simulationLogger.appendMethodCallStart("DpdSimulationTask.call", tmpIdDoInBackground);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation state PRE_PROCESSING">
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.PRE_PROCESSING);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Preparation of time step loop">
            if (!this.parameters.hasRestartInfo()) {
                this.initialiseVelocities();
            }
            
            ITimeStepCalculator tmpTimeStepCalculator = 
                this.factory.getTimeStepCalculator(
                    this.simulationOutput,
                    this.simulationLogger, 
                    this.parameters,
                    this.parallelizationInfo, 
                    this.randomNumberSeed,
                    this.particlePositionPool,
                    this.maximumNumberOfPositionCorrectionTrials
                );
            
            MoleculeRgValue[] tmpMoleculeRgValues = null;
            if (this.parameters.getChemicalSystemDescription().isRgCalculation()) {
                tmpMoleculeRgValues = new MoleculeRgValue[this.parameters.getChemicalSystemDescription().getRgCalculators().length];
                for (int i = 0; i < tmpMoleculeRgValues.length; i++) {
                    tmpMoleculeRgValues[i] = new MoleculeRgValue();
                }
            }
            
            IParticlePairInteractionCalculator tmpParticlePairNearestNeighborCalculator = null;
            if (this.parameters.getChemicalSystemDescription().isNearestNeighborParticleDetermination()) {
                // IMPORTANT: Use this.parameters.getChemicalSystemDescription().getNearestNeighborDistance() as cut-off length!
                tmpParticlePairNearestNeighborCalculator = 
                    this.factory.getParticlePairNearestNeighborCalculator(
                        this.simulationLogger, 
                        this.parameters.getChemicalSystemDescription().getBoxSize(), 
                        this.parameters.getSimulationDescription().getPeriodicBoundaries(), 
                        this.parameters.getChemicalSystemDescription().getNearestNeighborDistance(), 
                        this.parallelizationInfo, 
                        this.randomNumberSeed
                    );
            }
            HashMap<String, HashMap<String, Integer>> tmpBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap = null;
            HashMap<String, HashMap<String, Integer>> tmpBaseMoleculeParticleToNearestNeighborParticleFrequencyMap = null;
            HashMap<String, HashMap<String, Integer>> tmpBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap = null;
            HashMap<String, HashMap<String, Integer>> tmpBaseMoleculeToNearestNeighborMoleculeFrequencyMap = null;
            HashMap<String, HashMap<String, Integer>> tmpBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap = null;
            
            int tmpTimeStepNumber = this.parameters.getSimulationDescription().getTimeStepNumber();
            int tmpCurrentTimeStep = 0;
            if (this.parameters.hasRestartInfo()) {
                tmpTimeStepNumber = this.parameters.getRestartInfo().getAdditionalTimeStepNumber();
                tmpCurrentTimeStep = this.parameters.getRestartInfo().getLastTimeStep();
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation progress logging">
            this.simulationLogger.appendSimulationProgress("Preparations for time-step loop finished");
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation state TIME_STEP_INTEGRATION">
            long tmpStartCurrentTimeMillis = System.currentTimeMillis();
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.TIME_STEP_INTEGRATION);
            this.progressMonitor.setProgressInPercent(0);
            this.progressMonitor.setRemainingTime("---");
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Time-step loop">
            for (int i = 0; i < tmpTimeStepNumber; i++) {
                // <editor-fold defaultstate="collapsed" desc="Callable/thread interruption check">
                if (Thread.currentThread().isInterrupted()) {
                    this.isSimulationStopped = true;
                }
                // </editor-fold>
                ++tmpCurrentTimeStep;
                // <editor-fold defaultstate="collapsed" desc="Time step logging">
                long tmpIdTimeStep = this.simulationLogger.getId();
                this.simulationLogger.appendTimeStepStart(tmpCurrentTimeStep, tmpIdTimeStep);
                if (tmpCurrentTimeStep == 1 || tmpCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    this.simulationLogger.appendOutputTimeStepStart(tmpCurrentTimeStep, tmpIdTimeStep);
                }
                // </editor-fold>
                tmpTimeStepCalculator.calculate(tmpCurrentTimeStep);
                if (tmpCurrentTimeStep == 1 || tmpCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    // <editor-fold defaultstate="collapsed" desc="Property calculation">
                    if (this.parameters.getChemicalSystemDescription().isNearestNeighborParticleDetermination()) {
                        // <editor-fold defaultstate="collapsed" desc="Nearest-neighbor determination">
                        NearestNeighborManager tmpNearestNeighborManager = this.parameters.getChemicalSystemDescription().getNearestNeighborManager();
                        tmpNearestNeighborManager.reInitialize(
                            this.parameters.getParticleArrays().getNearestNeighborDistances(),
                            this.parameters.getParticleArrays().getNearestNeighborParticleIndices()
                        );
                        tmpParticlePairNearestNeighborCalculator.calculateParticlePairInteractions(
                            this.parameters.getParticleArrays().getR_x(),
                            this.parameters.getParticleArrays().getR_y(),
                            this.parameters.getParticleArrays().getR_z(),
                            this.parameters,
                            ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
                        );
                        tmpNearestNeighborManager.evaluateNearestNeighbors(
                            this.parameters.getParticleArrays().getNearestNeighborParticleIndices(),
                            this.parameters.getParticleArrays().getParticleTokens(),
                            this.parameters.getParticleArrays().getMoleculeNamesOfParticles(),
                            this.parameters.getParticleArrays().getMoleculeIndices()
                        );
                        tmpBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap = 
                            tmpNearestNeighborManager.getBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap();
                        tmpBaseMoleculeParticleToNearestNeighborParticleFrequencyMap = 
                            tmpNearestNeighborManager.getBaseMoleculeParticleToNearestNeighborParticleFrequencyMap();
                        tmpBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap = 
                            tmpNearestNeighborManager.getBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap();
                        tmpBaseMoleculeToNearestNeighborMoleculeFrequencyMap = 
                            tmpNearestNeighborManager.getBaseMoleculeToNearestNeighborMoleculeFrequencyMap();
                        tmpBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap = 
                            tmpNearestNeighborManager.getBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap();
                        // </editor-fold>
                    }
                    // Calculate potentials
                    tmpTimeStepCalculator.getPotentialAccumulator().accumulatePotentials(
                        this.parameters.getParticleArrays().getBondChunkArraysList(),
                        this.parameters.getParticleArrays().getR_x(),
                        this.parameters.getParticleArrays().getR_y(),
                        this.parameters.getParticleArrays().getR_z(),
                        this.parameters
                    );
                    // Upot
                    float tmpUpotDpd = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getDpdPotentialEnergyAdder().getSum();
                    float tmpUpotBonds = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getBondPotentialEnergyAdder().getSum();
                    float tmpUpotElectrostatics = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getElectrostaticsPotentialEnergyAdder().getSum();
                    float tmpUpotTotal = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getPotentialEnergyAdder().getSum();
                    // Surface tension
                    float tmpPressureX = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getPressureXAdder().getSum();
                    float tmpPressureY = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getPressureYAdder().getSum();
                    float tmpPressureZ = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getPressureZAdder().getSum();
                    float tmpDpdPressureX = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getDpdPressureXAdder().getSum();
                    float tmpDpdPressureY = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getDpdPressureYAdder().getSum();
                    float tmpDpdPressureZ = tmpTimeStepCalculator.getPotentialAccumulator().getExtendedAdderGroup().getDpdPressureZAdder().getSum();
                    float tmpBoxXLength = this.parameters.getChemicalSystemDescription().getBoxSize().getXLength();
                    float tmpBoxYLength = this.parameters.getChemicalSystemDescription().getBoxSize().getYLength();
                    float tmpBoxZLength = this.parameters.getChemicalSystemDescription().getBoxSize().getZLength();
                    float tmpSurfaceTensionAlongZ = (tmpPressureZ - 0.5f * (tmpPressureX + tmpPressureY)) / (tmpBoxXLength * tmpBoxYLength);
                    float tmpSurfaceTensionAlongX = (tmpPressureX - 0.5f * (tmpPressureZ + tmpPressureY)) / (tmpBoxZLength * tmpBoxYLength);
                    float tmpSurfaceTensionAlongY = (tmpPressureY - 0.5f * (tmpPressureX + tmpPressureZ)) / (tmpBoxXLength * tmpBoxZLength);
                    float tmpSurfaceTensionNorm = 
                        (float) FastMath.sqrt(tmpSurfaceTensionAlongX * tmpSurfaceTensionAlongX + tmpSurfaceTensionAlongY * tmpSurfaceTensionAlongY + tmpSurfaceTensionAlongZ * tmpSurfaceTensionAlongZ);
                    float tmpDpdSurfaceTensionAlongZ = (tmpDpdPressureZ - 0.5f * (tmpDpdPressureX + tmpDpdPressureY)) / (tmpBoxXLength * tmpBoxYLength);
                    float tmpDpdSurfaceTensionAlongX = (tmpDpdPressureX - 0.5f * (tmpDpdPressureZ + tmpDpdPressureY)) / (tmpBoxZLength * tmpBoxYLength);
                    float tmpDpdSurfaceTensionAlongY = (tmpDpdPressureY - 0.5f * (tmpDpdPressureX + tmpDpdPressureZ)) / (tmpBoxXLength * tmpBoxZLength);
                    float tmpDpdSurfaceTensionNorm = 
                        (float) FastMath.sqrt(tmpDpdSurfaceTensionAlongX * tmpDpdSurfaceTensionAlongX + tmpDpdSurfaceTensionAlongY * tmpDpdSurfaceTensionAlongY + tmpDpdSurfaceTensionAlongZ * tmpDpdSurfaceTensionAlongZ);
                    // Ukin
                    float tmpUkin = 
                        Utils.calculateUkin(
                            this.parameters.getParticleArrays().getV_x(),
                            this.parameters.getParticleArrays().getV_y(),
                            this.parameters.getParticleArrays().getV_z(),
                            this.parameters.getParticleArrays().getDpdMasses(),
                            this.parameters.getSimulationDescription().isDpdUnitMass()
                        );
                    // Utotal
                    float tmpUtotal = tmpUpotTotal + tmpUkin;
                    // Temperature
                    float tmpTemperature = tmpUkin / (ONE_A_HALF * (float) this.parameters.getParticleArrays().getV_x().length);
                    // Rg calculation if necessary
                    if (this.parameters.getChemicalSystemDescription().isRgCalculation()) {
                        RgCalculator[] tmpRgCalculators = this.parameters.getChemicalSystemDescription().getRgCalculators();
                        for (int k = 0; k < tmpRgCalculators.length; k++) {
                            float tmpMeanRadiusOfGyration = tmpRgCalculators[k].getMeanRadiusOfGyration(
                                this.parameters.getParticleArrays().getR_x(),
                                this.parameters.getParticleArrays().getR_y(),
                                this.parameters.getParticleArrays().getR_z(),
                                this.parameters.getParticleArrays().getMolarMasses()
                            );
                            tmpMoleculeRgValues[k].setMoleculeRgValue(tmpRgCalculators[k].getMoleculeName(), tmpMeanRadiusOfGyration);
                        }
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="Set output">
                    this.simulationOutput.setSimulationStepInformation(
                        tmpCurrentTimeStep, 
                        tmpTemperature, 
                        tmpUpotDpd,
                        tmpUpotBonds,
                        tmpUpotElectrostatics,
                        tmpUpotTotal, 
                        tmpUkin, 
                        tmpUtotal,
                        tmpSurfaceTensionAlongX,
                        tmpSurfaceTensionAlongY,
                        tmpSurfaceTensionAlongZ,
                        tmpSurfaceTensionNorm,
                        tmpDpdSurfaceTensionAlongX,
                        tmpDpdSurfaceTensionAlongY,
                        tmpDpdSurfaceTensionAlongZ,
                        tmpDpdSurfaceTensionNorm,
                        tmpMoleculeRgValues,
                        tmpBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap,
                        tmpBaseMoleculeParticleToNearestNeighborParticleFrequencyMap,
                        tmpBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap,
                        tmpBaseMoleculeToNearestNeighborMoleculeFrequencyMap,
                        tmpBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
                        Utils.getParticlePositions(this.parameters, this.particlePositionPool)
                    );
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Temperature              = " + String.valueOf(tmpTemperature));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Upot(DPD)                = " + String.valueOf(tmpUpotDpd));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Upot(Bonds)              = " + String.valueOf(tmpUpotBonds));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Upot(Electrostatics)     = " + String.valueOf(tmpUpotElectrostatics));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Upot(total)              = " + String.valueOf(tmpUpotTotal));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Ukin                     = " + String.valueOf(tmpUkin));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: U(total)                 = " + String.valueOf(tmpUtotal));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Surface Tension X        = " + String.valueOf(tmpSurfaceTensionAlongX));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Surface Tension Y        = " + String.valueOf(tmpSurfaceTensionAlongY));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Surface Tension Z        = " + String.valueOf(tmpSurfaceTensionAlongZ));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: Surface Tension Norm     = " + String.valueOf(tmpSurfaceTensionNorm));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: DPD Surface Tension X    = " + String.valueOf(tmpDpdSurfaceTensionAlongX));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: DPD Surface Tension Y    = " + String.valueOf(tmpDpdSurfaceTensionAlongY));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: DPD Surface Tension Z    = " + String.valueOf(tmpDpdSurfaceTensionAlongZ));
                    this.simulationLogger.appendIntermediateResults("DpdSimulationTask.call: DPD Surface Tension Norm = " + String.valueOf(tmpDpdSurfaceTensionNorm));
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="Particle force magnitude logging">
                    if (this.simulationLogger.isLogLevel(ILogLevel.PARTICLE)) {
                        tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().analyzeParticleForceMagnitudes(
                            this.parameters.getParticleArrays().getBondChunkArraysList(),
                            this.parameters.getParticleArrays().getR_x(),
                            this.parameters.getParticleArrays().getR_y(),
                            this.parameters.getParticleArrays().getR_z(),
                            this.parameters
                        );
                        float[] tmpMinMeanMaxDpdForceConservativeMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxDpdForceConservativeMagnitude();
                        if (tmpMinMeanMaxDpdForceConservativeMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: DPD conservat. F (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxDpdForceConservativeMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceConservativeMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceConservativeMagnitude[2])
                            );
                        }
                        float[] tmpMinMeanMaxDpdForceRandomMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxDpdForceRandomMagnitude();
                        if (tmpMinMeanMaxDpdForceRandomMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: DPD random F     (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxDpdForceRandomMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceRandomMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceRandomMagnitude[2])
                            );
                        }
                        float[] tmpMinMeanMaxDpdForceDissipativeMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxDpdForceDissipativeMagnitude();
                        if (tmpMinMeanMaxDpdForceDissipativeMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: DPD dissipat. F  (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[2])
                            );
                        }
                        float[] tmpMinMeanMaxBondForceMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxBondForceMagnitude();
                        if (tmpMinMeanMaxBondForceMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: Bond F           (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxBondForceMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxBondForceMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxBondForceMagnitude[2])
                            );
                        } else {
                            this.simulationLogger.appendParticleForceMagnitude("DpdSimulationTask.call: Bond F is not defined");
                        }
                        float[] tmpMinMeanMaxElectrostaticsForceMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxElectrostaticsForceMagnitude();
                        if (tmpMinMeanMaxElectrostaticsForceMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: Electrostatics F (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxElectrostaticsForceMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxElectrostaticsForceMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxElectrostaticsForceMagnitude[2])
                            );
                        } else {
                            this.simulationLogger.appendParticleForceMagnitude("DpdSimulationTask.call: Electrostatics F is not defined");
                        }
                        float[] tmpMinMeanMaxTotalForceMagnitude = tmpTimeStepCalculator.getParticleForceMagnitudeAccumulator().getMinMeanMaxTotalForceMagnitude();
                        if (tmpMinMeanMaxTotalForceMagnitude != null) {
                            this.simulationLogger.appendParticleForceMagnitude(
                                "DpdSimulationTask.call: Total F          (min/mean/max) = " + 
                                String.valueOf(tmpMinMeanMaxTotalForceMagnitude[0]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxTotalForceMagnitude[1]) +
                                " / " +
                                String.valueOf(tmpMinMeanMaxTotalForceMagnitude[2])
                            );
                        }
                        float[] tmpMinMeanMaxVelocityMagnitude =
                            Utils.calculateMinMeanMax(
                                this.parameters.getParticleArrays().getV_x(),
                                this.parameters.getParticleArrays().getV_y(),
                                this.parameters.getParticleArrays().getV_z()
                            );
                        this.simulationLogger.appendParticleForceMagnitude(
                            "DpdSimulationTask.call: Velocity v       (min/mean/max) = " + 
                            String.valueOf(tmpMinMeanMaxVelocityMagnitude[0]) +
                            " / " +
                            String.valueOf(tmpMinMeanMaxVelocityMagnitude[1]) +
                            " / " +
                            String.valueOf(tmpMinMeanMaxVelocityMagnitude[2])
                        );
                    }
                    // </editor-fold>
                }
                // <editor-fold defaultstate="collapsed" desc="Set progress monitor">
                float tmpSimulationProgressInPercent = ((float) i + 1.0f) / (float) tmpTimeStepNumber * HUNDRED;
                int tmpPercent = (int) tmpSimulationProgressInPercent;
                this.progressMonitor.setProgressInPercent(tmpPercent);
                long tmpProgressTimeMillis = System.currentTimeMillis() - tmpStartCurrentTimeMillis;
                long tmpTotalTimeMillis = (long) (HUNDRED / tmpSimulationProgressInPercent * (float) tmpProgressTimeMillis);
                this.progressMonitor.setRemainingTime(Utils.getFormattedTimePeriodString(tmpTotalTimeMillis - tmpProgressTimeMillis));
                // </editor-fold>
                // <editor-fold defaultstate="collapsed" desc="Time step logging and stop check">
                if (tmpCurrentTimeStep == 1 || tmpCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    this.simulationLogger.appendOutputTimeStepEnd(tmpCurrentTimeStep, tmpIdTimeStep);
                }
                this.simulationLogger.appendTimeStepEnd(tmpCurrentTimeStep, tmpIdTimeStep);
                if (tmpCurrentTimeStep == 1 || tmpCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    if (this.isSimulationStopped) {
                        this.progressMonitor.setSimulationStoppedFlag();
                        // <editor-fold defaultstate="collapsed" desc="Simulation progress logging">
                        this.simulationLogger.appendTimeStep("Time-step loop was manually stopped after tmpCurrentTimeStep = " + String.valueOf(tmpCurrentTimeStep));
                        // </editor-fold>
                        break;
                    }
                }
                // </editor-fold>
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation progress logging">
            this.simulationLogger.appendSimulationProgress("Time-step loop finished");
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation state POST_PROCESSING">
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.POST_PROCESSING);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Postprocessing after time step loop">
            RestartInfo tmpRestartInfo = 
                new RestartInfo(
                    tmpCurrentTimeStep,
                    this.parameters.getParticleArrays().getR_x(),
                    this.parameters.getParticleArrays().getR_y(),
                    this.parameters.getParticleArrays().getR_z(),
                    this.parameters.getParticleArrays().getV_x(),
                    this.parameters.getParticleArrays().getV_y(),
                    this.parameters.getParticleArrays().getV_z()
                );
            this.simulationOutput.setRestartInfo(tmpRestartInfo);
            tmpTimeStepCalculator.shutdownExecutorServices();
            if (tmpParticlePairNearestNeighborCalculator != null) {
                tmpParticlePairNearestNeighborCalculator.shutdownExecutorService();
            }
            if (!this.simulationOutput.finish()) {
                // <editor-fold defaultstate="collapsed" desc="Exception logging">
                this.simulationLogger.appendException("DpdSimulationTask.call", "WARNING: Simulation output did NOT finish successfully, output files may be corrupted.");
                // </editor-fold>
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulation progress logging">
            this.simulationLogger.appendSimulationProgress("Postprocessing for time-step loop finished");
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            this.simulationLogger.appendMethodCallEnd("DpdSimulationTask.call", tmpIdDoInBackground);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulator finished/stopped: Set progress monitor">
            this.progressMonitor.setProgressInPercent(100);
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.FINISHED_WITH_SUCCESS);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Finally release memory and finish logging">
            this.releaseMemoryAndFinishLogging();
            // </editor-fold>
            return true;
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("DpdSimulationTask.call", Utils.getStacktrace(anException));
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Simulator failed: Set progress monitor">
            this.progressMonitor.setSimulationState(IProgressMonitor.SimulationState.FINISHED_WITH_FAILURE);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Finally release memory and finish logging">
            this.releaseMemoryAndFinishLogging();
            // </editor-fold>
            return false;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    // <editor-fold defaultstate="collapsed" desc="- Initialisation related methods">
    /**
     * Sets simulation parameters
     * (No checks are performed)
     * 
     * @param aRestartInfo Restart info (may be null)
     * @param aSimulationInput Simulation input
     */
    private void setParameters(IInput aSimulationInput, RestartInfo aRestartInfo) {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("DpdSimulationTask.setParameters", tmpId);
        // </editor-fold>
        // IMPORTANT: Set parameters in the following order due to dependencies
        // <editor-fold defaultstate="collapsed" desc="1. SimulationCounts">
        SimulationCounts tmpSimulationCounts = new SimulationCounts(
            aSimulationInput.getParticleNumber()
        );
        this.simulationLogger.appendSimulationInit("SimulationCounts initialised");
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="2. ParticleTypes">
        ParticleTypes tmpParticleTypes = aSimulationInput.getParticleTypes();
        this.simulationLogger.appendSimulationInit("ParticleTypes initialised");
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="3. SimulationDescription">
        // Seed for random number generation:
        this.randomNumberSeed = new AtomicInteger(aSimulationInput.getRandomSeed());
        this.simulationLogger.appendSimulationInit("Random seed initialised");
        // Other simulation description parameters:
        SimulationDescription tmpSimulationDescription = new SimulationDescription(
            aSimulationInput.getTimeStepNumber(),
            aSimulationInput.getTimeStepLength(),
            aSimulationInput.getTimeStepFrequencyForOutput(),
            aSimulationInput.getInitialPotentialEnergyMinimizationStepNumber(),
            aSimulationInput.isInitialPotentialEnergyMinimizationWithAllForces(),
            aSimulationInput.isInitialPotentialEnergyMinimizationStepOutput(),
            aSimulationInput.getPeriodicBoundaries(),
            aSimulationInput.isDpdUnitMass(),
            aSimulationInput.getInitialVelocityScalingSteps()
        );
        this.simulationLogger.appendSimulationInit("SimulationDescription initialised");
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="4. InteractionDescription">
        float[][] tmpAij = aSimulationInput.getAij(tmpParticleTypes);
        // <editor-fold defaultstate="collapsed" desc="a(ij) logging">
        if (this.simulationLogger.isLogLevel(ILogLevel.A_IJ)) {
            for (int i = 0; i < tmpParticleTypes.getParticleTypeNumber(); i++) {
                for (int k = 0; k < tmpParticleTypes.getParticleTypeNumber(); k++) {
                    this.simulationLogger.appendAij("a(" + tmpParticleTypes.getParticleToken(i) + ", " + tmpParticleTypes.getParticleToken(k) + ") = " + String.valueOf(tmpAij[i][k]));
                }
            }
        }
        // </editor-fold>
        Electrostatics tmpElectrostatics = aSimulationInput.getElectrostatics();
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging for electrostatics">
        if (tmpElectrostatics != null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Electrostatics parameters defined");
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Electrostatics charge distribution type = " + tmpElectrostatics.getChargeDistributionType().toString());
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Electrostatics splitting type           = " + tmpElectrostatics.getSplittingType().toString());
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No electrostatics parameters");
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging for graviational acceleration">
        if (aSimulationInput.getGravitationalAcceleration().isGravitationalAcceleration()) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Gravitational acceleration is defined");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No gravitational acceleration ");
        }
        // </editor-fold>
        InteractionDescription tmpInteractionDescription = new InteractionDescription(
            aSimulationInput.getTemperature(),
            aSimulationInput.getDpdSigma(),
            aSimulationInput.isGaussianRandomDpdForce(),
            tmpParticleTypes.getParticleTypeNumber(),
            tmpAij,
            tmpElectrostatics,
            aSimulationInput.getGravitationalAcceleration(),
            tmpSimulationDescription
        );
        this.simulationLogger.appendSimulationInit("InteractionDescription initialised");
        this.simulationLogger.appendIntermediateResults("DPD Temperature          = " + String.valueOf(tmpInteractionDescription.getTemperature()));
        this.simulationLogger.appendIntermediateResults("DPD Sigma                = " + String.valueOf(tmpInteractionDescription.getDpdSigma()));
        this.simulationLogger.appendIntermediateResults("DPD Sigma/SQRT(TimeStep) = " + String.valueOf(tmpInteractionDescription.getDpdSigmaDivRootTimeStepLength()));
        this.simulationLogger.appendIntermediateResults("DPD Gamma                = " + String.valueOf(tmpInteractionDescription.getDpdGamma()));
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="5. ParticleArrays and ChemicalSystemDescription">
        MoleculeDescription[] tmpMoleculeDescriptions = aSimulationInput.getMoleculeDescriptions();
        // Initialise MoleculeTypes arrays
        HashMap<String, Integer> tmpMoleculeNameToIndexMap = new HashMap<>(tmpMoleculeDescriptions.length);
        String[] tmpMoleculeNames = new String[tmpMoleculeDescriptions.length];
        int[] tmpTotalMoleculeParticleNumbers = new int[tmpMoleculeDescriptions.length];
        int[] tmpSingleMoleculeParticleNumbers = new int[tmpMoleculeDescriptions.length];
        int[] tmpFirstIndices = new int[tmpMoleculeDescriptions.length];
        int[] tmpLastIndices = new int[tmpMoleculeDescriptions.length];
        // Initialise particle arrays
        float[] tmpR_x = new float[tmpSimulationCounts.getParticleNumber()];
        float[] tmpR_y = new float[tmpSimulationCounts.getParticleNumber()];
        float[] tmpR_z = new float[tmpSimulationCounts.getParticleNumber()];
        String[] tmpParticleTokens = new String[tmpSimulationCounts.getParticleNumber()];
        String[] tmpMoleculeNamesOfParticles = new String[tmpSimulationCounts.getParticleNumber()];
        int[] tmpParticleTypeIndices = new int[tmpSimulationCounts.getParticleNumber()];
        int[] tmpMoleculeTypeIndices = new int[tmpSimulationCounts.getParticleNumber()];
        int[] tmpMoleculeIndices = new int[tmpSimulationCounts.getParticleNumber()];
        float[] tmpParticleCharges = new float[tmpSimulationCounts.getParticleNumber()];
        float[] tmpParticleDpdMasses = new float[tmpSimulationCounts.getParticleNumber()];
        float[] tmpParticleMolarMasses = new float[tmpSimulationCounts.getParticleNumber()];
        LinkedList<Integer> tmpChargedParticleIndexList = new LinkedList<>();
        // Check bond existence
        HarmonicBondChunkGenerator tmpBondChunkGenerator = null;
        HashMap<String, ParticlePairHarmonicBond> tmpParticlePairBondMap = null;
        for (MoleculeDescription tmpMoleculeDescription : tmpMoleculeDescriptions) {
            if (tmpMoleculeDescription.getSingleMoleculeParticleNumber() > 1) {
                tmpBondChunkGenerator = new HarmonicBondChunkGenerator();
                tmpParticlePairBondMap = aSimulationInput.getParticlePairBondMap();
                break;
            }
        }
        // Loop over all molecule descriptions
        int tmpMoleculeIndex = 0;
        int tmpFirstIndex = 0;
        for (int i = 0; i < tmpMoleculeDescriptions.length; i++) {
            MoleculeDescription tmpMoleculeDescription = tmpMoleculeDescriptions[i];

            // tmpFirstIndex = 0, getTotalMoleculeParticleNumber() = 3 -> tmpLastIndex = 2
            int tmpLastIndex = tmpFirstIndex + tmpMoleculeDescription.getTotalMoleculeParticleNumber() - 1;

            String tmpMoleculeName = tmpMoleculeDescription.getMoleculeName();
            tmpMoleculeNameToIndexMap.put(tmpMoleculeName, i);
            tmpMoleculeNames[i] = tmpMoleculeName;
            tmpTotalMoleculeParticleNumbers[i] = tmpMoleculeDescription.getTotalMoleculeParticleNumber();
            tmpSingleMoleculeParticleNumbers[i] = tmpMoleculeDescription.getSingleMoleculeParticleNumber();
            tmpFirstIndices[i] = tmpFirstIndex;
            tmpLastIndices[i] = tmpLastIndex;

            HashMap<Integer, Integer> tmpBackBoneIndexToArrayIndexMap = new HashMap<>(tmpMoleculeDescription.getParticleBackboneIndices().length);
            // IMPORTANT: Set to 1 since increment is at the end of the loop
            int tmpMoleculeParticleCounter = 1;
            
            int tmpCurrentIndex = tmpFirstIndex;
            for (int j = 0; j < tmpMoleculeDescription.getTotalMoleculeParticleNumber(); j++) {
                tmpR_x[tmpCurrentIndex] = tmpMoleculeDescription.getR_x()[j];
                tmpR_y[tmpCurrentIndex] = tmpMoleculeDescription.getR_y()[j];
                tmpR_z[tmpCurrentIndex] = tmpMoleculeDescription.getR_z()[j];

                if (tmpMoleculeDescription.getParticleBackboneIndices()[j] > 0) {
                    tmpBackBoneIndexToArrayIndexMap.put(tmpMoleculeDescription.getParticleBackboneIndices()[j], tmpCurrentIndex);
                }
                
                String tmpParticleToken = tmpMoleculeDescription.getParticleTokens()[j];
                tmpParticleTokens[tmpCurrentIndex] = tmpParticleToken;
                tmpMoleculeNamesOfParticles[tmpCurrentIndex] = tmpMoleculeName;
                int tmpParticleTypeIndex = tmpParticleTypes.getIndex(tmpParticleToken);
                tmpParticleTypeIndices[tmpCurrentIndex] = tmpParticleTypeIndex;
                tmpMoleculeTypeIndices[tmpCurrentIndex] = i;
                tmpParticleCharges[tmpCurrentIndex] = tmpParticleTypes.getCharge(tmpParticleTypeIndex);
                if (tmpParticleTypes.getCharge(tmpParticleTypeIndex) != 0) {
                    tmpChargedParticleIndexList.add(tmpCurrentIndex);
                }
                if (tmpSimulationDescription.isDpdUnitMass()) {
                    tmpParticleDpdMasses[tmpCurrentIndex] = 1.0f;
                } else {
                    tmpParticleDpdMasses[tmpCurrentIndex] = tmpParticleTypes.getMolarMass(tmpParticleTypeIndex)/tmpParticleTypes.getMinimumMolarMass();
                }
                tmpParticleMolarMasses[tmpCurrentIndex] = tmpParticleTypes.getMolarMass(tmpParticleTypeIndex);
                
                // Set bonds
                int[] tmpSingleParticleBondOffsets = tmpMoleculeDescription.getBondOffsets()[j];
                if (tmpSingleParticleBondOffsets != null) {
                    for (int k = 0; k < tmpSingleParticleBondOffsets.length; k++) {
                        // Only backward bonds are evaluated to avoid float counting
                        if (tmpSingleParticleBondOffsets[k] < 0) {
                            String tmpParticleToken1 = tmpParticleToken;
                            String tmpParticleToken2 = tmpParticleTokens[tmpCurrentIndex + tmpSingleParticleBondOffsets[k]];
                            String tmpParticleTokenPairKey = Utils.getParticleTokenPairKey(tmpParticleToken1, tmpParticleToken2);
                            ParticlePairHarmonicBond tmpParticlePairBond = tmpParticlePairBondMap.get(tmpParticleTokenPairKey);
                            tmpBondChunkGenerator.addBond(
                                new HarmonicBond(
                                    tmpCurrentIndex,
                                    tmpCurrentIndex + tmpSingleParticleBondOffsets[k],
                                    tmpParticlePairBond.getBondLength(),
                                    tmpParticlePairBond.getForceConstant(),
                                    tmpParticlePairBond.getHarmonicBondBehaviour()
                                )
                            );
                        }
                    }
                }
                
                // Set backbone bonds of single molecule
                if (tmpMoleculeDescription.getBackboneBondNumber() > 0 && tmpMoleculeParticleCounter == tmpMoleculeDescription.getSingleMoleculeParticleNumber()) {
                    for (HarmonicBond tmpSingleBackboneBond : tmpMoleculeDescription.getBackboneBonds()) {
                        tmpBondChunkGenerator.addBond(
                            new HarmonicBond(
                                tmpBackBoneIndexToArrayIndexMap.get(tmpSingleBackboneBond.getIndex1()),
                                tmpBackBoneIndexToArrayIndexMap.get(tmpSingleBackboneBond.getIndex2()),
                                tmpSingleBackboneBond.getBondLength(),
                                tmpSingleBackboneBond.getForceConstant(),
                                tmpSingleBackboneBond.getHarmonicBondBehaviour()
                            )
                        );
                    }
                    // NOTE: The clear method is not necessary (since old entries are simply replaced by new ones with put-operation) 
                    //       but better programming style
                    tmpBackBoneIndexToArrayIndexMap.clear();
                }
                
                // Set molecule index
                tmpMoleculeIndices[tmpCurrentIndex] = tmpMoleculeIndex;
                if (tmpMoleculeParticleCounter == tmpMoleculeDescription.getSingleMoleculeParticleNumber()) {
                    // IMPORTANT: Set to 0 since increment is at the end of the loop
                    tmpMoleculeParticleCounter = 0;
                    tmpMoleculeIndex++;
                }                
                
                tmpMoleculeParticleCounter++;
                tmpCurrentIndex++;
            }
            tmpFirstIndex = tmpLastIndex + 1;
        }
        LinkedList<HarmonicBondChunkArrays> tmpBondChunkArraysList = null;
        if (tmpBondChunkGenerator != null && tmpBondChunkGenerator.hasBonds()) {
            tmpBondChunkArraysList = tmpBondChunkGenerator.getBondChunkArraysList();
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpBondChunkGenerator != null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Bond number                    = " + String.valueOf(tmpBondChunkGenerator.getBondNumber()));
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Bond chunk number              = " + String.valueOf(tmpBondChunkGenerator.getBondChunkNumber()));
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Bond info                      : " + tmpBondChunkGenerator.getBondInfo());
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No bonds");
        }
        // </editor-fold>
        // Set charged particles
        int[] tmpChargedParticleIndices = null;
        if (!tmpChargedParticleIndexList.isEmpty()) {
            tmpChargedParticleIndices = new int[tmpChargedParticleIndexList.size()];
            int tmpIndex = 0;
            for (int tmpChargedParticleIndex : tmpChargedParticleIndexList) {
                tmpChargedParticleIndices[tmpIndex++] = tmpChargedParticleIndex;
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpChargedParticleIndices != null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Charged particle number        = " + String.valueOf(tmpChargedParticleIndices.length));
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No charged particles");
        }
        // </editor-fold>
        // Calculated molar masses of molecules
        float[] tmpMoleculeMolarMasses = new float[tmpMoleculeNames.length];
        for (int i = 0; i < tmpMoleculeNames.length; i++) {
            float tmpMoleculeMolarMass = 0;
            for (int j = 0; j < tmpSingleMoleculeParticleNumbers[i]; j++) {
                tmpMoleculeMolarMass += tmpParticleMolarMasses[tmpFirstIndices[i] + j];
            }
            tmpMoleculeMolarMasses[i] = tmpMoleculeMolarMass;
        }
        // Set molecule types
        MoleculeTypes tmpMoleculeTypes = 
            new MoleculeTypes(
                tmpMoleculeNameToIndexMap, 
                tmpMoleculeNames, 
                tmpMoleculeMolarMasses,
                tmpTotalMoleculeParticleNumbers, 
                tmpSingleMoleculeParticleNumbers, 
                tmpFirstIndices, 
                tmpLastIndices
            );
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Molecule number                = " + String.valueOf(tmpMoleculeNames.length));
        // </editor-fold>
        // Set molecule fixations if necessary
        MoleculeFixationDescription[] tmpMoleculeFixationDescriptions = aSimulationInput.getMoleculeFixationDescriptions();
        MoleculeFixationInfo[] tmpMoleculeFixationInfos = null;
        if (tmpMoleculeFixationDescriptions != null) {
            tmpMoleculeFixationInfos = new MoleculeFixationInfo[tmpMoleculeFixationDescriptions.length];
            for (int i = 0; i < tmpMoleculeFixationDescriptions.length; i++) {
                tmpMoleculeFixationInfos[i] =
                    new MoleculeFixationInfo(
                        tmpMoleculeFixationDescriptions[i].getMoleculeName(),
                        tmpMoleculeFixationDescriptions[i].isFixedX(),
                        tmpMoleculeFixationDescriptions[i].isFixedY(),
                        tmpMoleculeFixationDescriptions[i].isFixedZ(),
                        tmpMoleculeFixationDescriptions[i].getMaxTimeStep(),
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeFixationDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getLastIndex(tmpMoleculeFixationDescriptions[i].getMoleculeName())
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpMoleculeFixationInfos == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule fixations");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Molecule fixations             = " + String.valueOf(tmpMoleculeFixationInfos.length));
        }
        // </editor-fold>
        // Set molecule boundaries if necessary
        MoleculeBoundaryDescription[] tmpMoleculeBoundaryDescriptions = aSimulationInput.getMoleculeBoundaryDescriptions();
        MoleculeBoundaryInfo[] tmpMoleculeBoundaryInfos = null;
        if (tmpMoleculeBoundaryDescriptions != null) {
            tmpMoleculeBoundaryInfos = new MoleculeBoundaryInfo[tmpMoleculeBoundaryDescriptions.length];
            for (int i = 0; i < tmpMoleculeBoundaryDescriptions.length; i++) {
                tmpMoleculeBoundaryInfos[i] =
                    new MoleculeBoundaryInfo(
                        tmpMoleculeBoundaryDescriptions[i].getMoleculeName(),
                        tmpMoleculeBoundaryDescriptions[i].isBoundaryX(),
                        tmpMoleculeBoundaryDescriptions[i].getXmin(),
                        tmpMoleculeBoundaryDescriptions[i].getXmax(),
                        tmpMoleculeBoundaryDescriptions[i].isBoundaryY(),
                        tmpMoleculeBoundaryDescriptions[i].getYmin(),
                        tmpMoleculeBoundaryDescriptions[i].getYmax(),
                        tmpMoleculeBoundaryDescriptions[i].isBoundaryZ(),
                        tmpMoleculeBoundaryDescriptions[i].getZmin(),
                        tmpMoleculeBoundaryDescriptions[i].getZmax(),
                        tmpMoleculeBoundaryDescriptions[i].getMaxTimeStep(),
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeBoundaryDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getLastIndex(tmpMoleculeBoundaryDescriptions[i].getMoleculeName())
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpMoleculeBoundaryInfos == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule boundaries");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Molecule boundaries            = " + String.valueOf(tmpMoleculeBoundaryInfos.length));
        }
        // </editor-fold>
        // Set molecule spheres if necessary
        MoleculeSphereDescription[] tmpMoleculeSphereDescriptions = aSimulationInput.getMoleculeSphereDescriptions();
        MoleculeSphereInfo[] tmpMoleculeSphereInfos = null;
        if (tmpMoleculeSphereDescriptions != null) {
            tmpMoleculeSphereInfos = new MoleculeSphereInfo[tmpMoleculeSphereDescriptions.length];
            for (int i = 0; i < tmpMoleculeSphereDescriptions.length; i++) {
                tmpMoleculeSphereInfos[i] =
                    new MoleculeSphereInfo(
                        tmpMoleculeSphereDescriptions[i].getMoleculeName(),
                        tmpMoleculeSphereDescriptions[i].isExclusiveSphere(),
                        tmpMoleculeSphereDescriptions[i].getSphereCenterX(),
                        tmpMoleculeSphereDescriptions[i].getSphereCenterY(),
                        tmpMoleculeSphereDescriptions[i].getSphereCenterZ(),
                        tmpMoleculeSphereDescriptions[i].getSphereRadius(),
                        tmpMoleculeSphereDescriptions[i].getMaxTimeStep(),
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeSphereDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getLastIndex(tmpMoleculeSphereDescriptions[i].getMoleculeName())
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpMoleculeSphereInfos == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule exlusion/inclusion spheres");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Molecule incl./excl. spheres   = " + String.valueOf(tmpMoleculeSphereInfos.length));
        }
        // </editor-fold>
        // Set molecule velocity fixations if necessary
        MoleculeVelocityFixationDescription[] tmpMoleculeVelocityFixationDescriptions = aSimulationInput.getMoleculeVelocityFixationDescriptions();
        MoleculeVelocityFixationInfo[] tmpMoleculeVelocityFixationInfos = null;
        if (tmpMoleculeVelocityFixationDescriptions != null) {
            tmpMoleculeVelocityFixationInfos = new MoleculeVelocityFixationInfo[tmpMoleculeVelocityFixationDescriptions.length];
            for (int i = 0; i < tmpMoleculeVelocityFixationDescriptions.length; i++) {
                tmpMoleculeVelocityFixationInfos[i] =
                    new MoleculeVelocityFixationInfo(
                        tmpMoleculeVelocityFixationDescriptions[i].getMoleculeName(),
                        tmpMoleculeVelocityFixationDescriptions[i].isFixedX(),
                        tmpMoleculeVelocityFixationDescriptions[i].isFixedY(),
                        tmpMoleculeVelocityFixationDescriptions[i].isFixedZ(),
                        tmpMoleculeVelocityFixationDescriptions[i].getVelocityX(),
                        tmpMoleculeVelocityFixationDescriptions[i].getVelocityY(),
                        tmpMoleculeVelocityFixationDescriptions[i].getVelocityZ(),
                        tmpMoleculeVelocityFixationDescriptions[i].getMaxTimeStep(),
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeVelocityFixationDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getLastIndex(tmpMoleculeVelocityFixationDescriptions[i].getMoleculeName())
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpMoleculeVelocityFixationInfos == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule velocity fixations");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Velocity fixations             = " + String.valueOf(tmpMoleculeVelocityFixationInfos.length));
        }
        // </editor-fold>
        // Set molecule accelerations if necessary
        MoleculeAccelerationDescription[] tmpMoleculeAccelerationDescriptions = aSimulationInput.getMoleculeAccelerationDescriptions();
        MoleculeAccelerationInfo[] tmpMoleculeAccelerationInfos = null;
        if (tmpMoleculeAccelerationDescriptions != null) {
            tmpMoleculeAccelerationInfos = new MoleculeAccelerationInfo[tmpMoleculeAccelerationDescriptions.length];
            for (int i = 0; i < tmpMoleculeAccelerationDescriptions.length; i++) {
                tmpMoleculeAccelerationInfos[i] =
                    new MoleculeAccelerationInfo(
                        tmpMoleculeAccelerationDescriptions[i].getMoleculeName(),
                        tmpMoleculeAccelerationDescriptions[i].getAccelerationX(),
                        tmpMoleculeAccelerationDescriptions[i].getAccelerationY(),
                        tmpMoleculeAccelerationDescriptions[i].getAccelerationZ(),
                        tmpMoleculeAccelerationDescriptions[i].getFrequency(),
                        tmpMoleculeAccelerationDescriptions[i].getMaxTimeStep(),
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeAccelerationDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getLastIndex(tmpMoleculeAccelerationDescriptions[i].getMoleculeName())
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpMoleculeAccelerationInfos == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule accelerations");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Velocity accelerations         = " + String.valueOf(tmpMoleculeAccelerationInfos.length));
        }
        // </editor-fold>
        // Box size 
        BoxSize tmpBoxSize = aSimulationInput.getBoxSize();
        // Set molecule particle radius of gyration (Rg) calculation if necessary
        MoleculeRgCalculationDescription[] tmpMoleculeRgCalculationDescriptions = aSimulationInput.getMoleculeRgCalculationDescriptions();
        RgCalculator[] tmpRgCalculators = null;
        if (tmpMoleculeRgCalculationDescriptions != null) {
            tmpRgCalculators = new RgCalculator[tmpMoleculeRgCalculationDescriptions.length];
            for (int i = 0; i < tmpMoleculeRgCalculationDescriptions.length; i++) {
                tmpRgCalculators[i] =
                    new RgCalculator(
                        tmpMoleculeTypes.getFirstIndex(tmpMoleculeRgCalculationDescriptions[i].getMoleculeName()),
                        tmpMoleculeRgCalculationDescriptions[i].getMoleculeName(),
                        tmpMoleculeTypes.getMoleculeMolarMass(tmpMoleculeRgCalculationDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getSingleMoleculeParticleNumber(tmpMoleculeRgCalculationDescriptions[i].getMoleculeName()),
                        tmpMoleculeTypes.getTotalMoleculeParticleNumber(tmpMoleculeRgCalculationDescriptions[i].getMoleculeName()),
                        tmpBoxSize,
                        tmpSimulationDescription.getPeriodicBoundaries()
                    );
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpRgCalculators == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No molecule Rg calculations");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: Rg calculations                = " + String.valueOf(tmpRgCalculators.length));
        }
        // </editor-fold>
        // Set nearest-neigbor particle manager if necessary
        NearestNeighborBaseParticleDescription[] tmpNearestNeighborParticleDescriptions = aSimulationInput.getNearestNeighborParticleDescriptions();
        NearestNeighborManager tmpNearestNeighborManager = null;
        boolean[] tmpIsNearestNeighborBaseParticleDeterminations = null;
        float[] tmpNearestNeighborDistances = null;
        int[] tmpNearestNeighborParticleIndices = null;
        if (tmpNearestNeighborParticleDescriptions != null) {
            // Initialize nearest-neighbor related arrays
            tmpIsNearestNeighborBaseParticleDeterminations = new boolean[tmpSimulationCounts.getParticleNumber()];
            Arrays.fill(tmpIsNearestNeighborBaseParticleDeterminations, false);
            tmpNearestNeighborDistances = new float[tmpSimulationCounts.getParticleNumber()];
            tmpNearestNeighborParticleIndices = new int[tmpSimulationCounts.getParticleNumber()];
            // Initialize nearest-neighbor manager
            tmpNearestNeighborManager = new NearestNeighborManager();
            for (NearestNeighborBaseParticleDescription tmpNearestNeighborParticleDescription : tmpNearestNeighborParticleDescriptions) {
                String tmpMoleculeName = tmpNearestNeighborParticleDescription.getMoleculeName();
                tmpFirstIndex = tmpMoleculeTypes.getFirstIndex(tmpMoleculeName);
                int tmpLastIndex =  tmpMoleculeTypes.getLastIndex(tmpMoleculeName);
                String tmpParticleToken = tmpNearestNeighborParticleDescription.getBaseParticleToken();
                LinkedList<Integer> tmpParticleIndexList = new LinkedList<>();
                for (int i = tmpFirstIndex; i <= tmpLastIndex; i++) {
                    if (tmpParticleTokens[i].equals(tmpParticleToken)) {
                        tmpParticleIndexList.add(i);
                        tmpIsNearestNeighborBaseParticleDeterminations[i] = true;
                    }
                }
                tmpNearestNeighborParticleDescription.setBaseParticleIndexList(tmpParticleIndexList);
                tmpNearestNeighborManager.addNearestNeighborParticleDescription(tmpNearestNeighborParticleDescription);
            }
            tmpNearestNeighborManager.consolidate();
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (tmpNearestNeighborManager == null) {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: No nearest-neigbor determination");
        } else {
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.setParameters: " + String.valueOf(tmpNearestNeighborParticleDescriptions.length) + " nearest-neigbor determinations defined");
        }
        // </editor-fold>
        // Set ParticleArrays instance
        ParticleArrays tmpParticleArrays = null;
        if (aRestartInfo == null) {
            tmpParticleArrays = 
                new ParticleArrays(
                    tmpR_x,
                    tmpR_y,
                    tmpR_z,
                    new float[tmpSimulationCounts.getParticleNumber()],
                    new float[tmpSimulationCounts.getParticleNumber()],
                    new float[tmpSimulationCounts.getParticleNumber()],
                    tmpParticleTokens,
                    tmpMoleculeNamesOfParticles,
                    tmpParticleTypeIndices,
                    tmpMoleculeTypeIndices,
                    tmpMoleculeIndices,                        
                    tmpParticleCharges, 
                    tmpParticleDpdMasses,
                    tmpParticleMolarMasses,
                    tmpBondChunkArraysList,
                    tmpChargedParticleIndices,
                    tmpIsNearestNeighborBaseParticleDeterminations,
                    tmpNearestNeighborDistances,
                    tmpNearestNeighborParticleIndices
                );
            this.simulationLogger.appendSimulationInit("ParticleArrays initialised WITHOUT aRestartInfo");
        } else {
            tmpParticleArrays = 
                new ParticleArrays(
                    aRestartInfo.getR_x(),
                    aRestartInfo.getR_y(),
                    aRestartInfo.getR_z(),
                    aRestartInfo.getV_x(),
                    aRestartInfo.getV_y(),
                    aRestartInfo.getV_z(),
                    tmpParticleTokens, 
                    tmpMoleculeNamesOfParticles,
                    tmpParticleTypeIndices,
                    tmpMoleculeTypeIndices,
                    tmpMoleculeIndices,                        
                    tmpParticleCharges, 
                    tmpParticleDpdMasses,
                    tmpParticleMolarMasses,
                    tmpBondChunkArraysList,
                    tmpChargedParticleIndices,
                    tmpIsNearestNeighborBaseParticleDeterminations,
                    tmpNearestNeighborDistances,
                    tmpNearestNeighborParticleIndices
                );
            this.simulationLogger.appendSimulationInit("ParticleArrays initialised WITH aRestartInfo");
        }
        // Set chemical system description
        ChemicalSystemDescription tmpChemicalSystemDescription =
            new ChemicalSystemDescription(
                tmpMoleculeTypes,
                tmpBoxSize,
                tmpMoleculeFixationInfos,
                tmpMoleculeBoundaryInfos,
                tmpMoleculeSphereInfos,
                tmpMoleculeVelocityFixationInfos,
                tmpMoleculeAccelerationInfos,
                tmpRgCalculators,
                tmpNearestNeighborManager,
                aSimulationInput.getNearestNeighborDistance()
            );
        this.simulationLogger.appendSimulationInit("ChemicalSystemDescription initialised");
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="6. Test objects">
        // Test objects are NOT defined in operational DpdSimulationTask
        TestObjects tmpTestObjects = 
            new TestObjects(
                null,
                null,
                null,
                null
            );
        // </editor-fold>
        this.parameters = new Parameters(
            aRestartInfo, 
            tmpParticleArrays, 
            tmpParticleTypes, 
            tmpChemicalSystemDescription, 
            tmpInteractionDescription, 
            tmpSimulationDescription, 
            tmpSimulationCounts, 
            tmpTestObjects
        );
        this.simulationLogger.appendSimulationInit("Parameters initialised");
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("DpdSimulationTask.setParameters", tmpId);
        // </editor-fold>
    }
    
    /**
     * NO restart: Initialises velocities 
     * this.parameters.getParticleArrays().getV_x(), 
     * this.parameters.getParticleArrays().getV_y() and 
     * this.parameters.getParticleArrays().getV_z()
     */
    private void initialiseVelocities() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("DpdSimulationTask.initialiseVelocities (NO restart)", tmpId);
        // </editor-fold>
        IRandom tmpRandomX = this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet());
        IRandom tmpRandomY = this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet());
        IRandom tmpRandomZ = this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet());
        for (int i = 0; i < this.parameters.getParticleArrays().getV_x().length; i++) {
            this.parameters.getParticleArrays().getV_x()[i] = tmpRandomX.nextGaussian();
            this.parameters.getParticleArrays().getV_y()[i] = tmpRandomY.nextGaussian();
            this.parameters.getParticleArrays().getV_z()[i] = tmpRandomZ.nextGaussian();
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        if (this.simulationLogger.isLogLevel(ILogLevel.QUANTITY)) {
            float tmpUkin = 
                Utils.calculateUkin(
                    this.parameters.getParticleArrays().getV_x(),
                    this.parameters.getParticleArrays().getV_y(),
                    this.parameters.getParticleArrays().getV_z(),
                    this.parameters.getParticleArrays().getDpdMasses(),
                    this.parameters.getSimulationDescription().isDpdUnitMass()
                );
            this.simulationLogger.appendIntermediateResults("DpdSimulationTask.initialiseVelocities, tmpUkin = " + String.valueOf(tmpUkin));
        }
        // </editor-fold>
        float tmpVelocityScaleFactor = Utils.scale_v(
            this.parameters.getParticleArrays().getV_x(),
            this.parameters.getParticleArrays().getV_y(),
            this.parameters.getParticleArrays().getV_z(),
            this.parameters.getParticleArrays().getDpdMasses(),
            this.parameters.getInteractionDescription().getTemperature(),
            this.parameters.getSimulationDescription().isDpdUnitMass()
        );
        // <editor-fold defaultstate="collapsed" desc="Velocity scale factor logging">
        this.simulationLogger.appendVelocityScaleFactor("DpdSimulationTask.initialiseVelocities, velocity scale factor = " + String.valueOf(tmpVelocityScaleFactor));
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("DpdSimulationTask.initialiseVelocities (NO restart)", tmpId);
        // </editor-fold>
    }

    /**
     * Corrects molecule boundaries if necessary.
     */
    private void correctMoleculeBoundaries() {
        if (this.parameters.getChemicalSystemDescription().getMoleculeBoundaryInfos() != null) {
            float[] tmpR_x = this.parameters.getParticleArrays().getR_x();
            float[] tmpR_y = this.parameters.getParticleArrays().getR_y();
            float[] tmpR_z = this.parameters.getParticleArrays().getR_z();
            for (MoleculeBoundaryInfo tmpMoleculeBoundaryInfo: this.parameters.getChemicalSystemDescription().getMoleculeBoundaryInfos()) {
                if (tmpMoleculeBoundaryInfo.isBoundaryX()) {
                    float tmpXmin = tmpMoleculeBoundaryInfo.getXmin();
                    float tmpXmax = tmpMoleculeBoundaryInfo.getXmax();
                    boolean tmpIsMinSmallerMax = tmpXmin < tmpXmax;
                    boolean tmpIsMinGreaterMax = tmpXmin > tmpXmax;
                    for (int i = tmpMoleculeBoundaryInfo.getFirstIndex(); i < tmpMoleculeBoundaryInfo.getExclusiveLastIndex(); i++) {
                        if (tmpIsMinSmallerMax) {
                            if (tmpR_x[i] < tmpXmin) tmpR_x[i] = tmpXmin;
                            if (tmpR_x[i] > tmpXmax) tmpR_x[i] = tmpXmax;
                        } else if (tmpIsMinGreaterMax) {
                            if (tmpR_x[i] < tmpXmin && tmpR_x[i] > tmpXmax) {
                                if (tmpXmin - tmpR_x[i] <= tmpR_x[i] - tmpXmax) tmpR_x[i] = tmpXmin;
                                if (tmpXmin - tmpR_x[i] > tmpR_x[i] - tmpXmax) tmpR_x[i] = tmpXmax;
                            }
                        }
                        // Min = Max: Do nothing!
                    }
                }
                if (tmpMoleculeBoundaryInfo.isBoundaryY()) {
                    float tmpYmin = tmpMoleculeBoundaryInfo.getYmin();
                    float tmpYmax = tmpMoleculeBoundaryInfo.getYmax();
                    boolean tmpIsMinSmallerMax = tmpYmin < tmpYmax;
                    boolean tmpIsMinGreaterMax = tmpYmin > tmpYmax;
                    for (int i = tmpMoleculeBoundaryInfo.getFirstIndex(); i < tmpMoleculeBoundaryInfo.getExclusiveLastIndex(); i++) {
                        if (tmpIsMinSmallerMax) {
                            if (tmpR_y[i] < tmpYmin) tmpR_y[i] = tmpYmin;
                            if (tmpR_y[i] > tmpYmax) tmpR_y[i] = tmpYmax;
                        } else if (tmpIsMinGreaterMax) {
                            if (tmpR_y[i] < tmpYmin && tmpR_y[i] > tmpYmax) {
                                if (tmpYmin - tmpR_y[i] <= tmpR_y[i] - tmpYmax) tmpR_y[i] = tmpYmin;
                                if (tmpYmin - tmpR_y[i] > tmpR_y[i] - tmpYmax) tmpR_y[i] = tmpYmax;
                            }
                        }
                        // Min = Max: Do nothing!
                    }
                }
                if (tmpMoleculeBoundaryInfo.isBoundaryZ()) {
                    float tmpZmin = tmpMoleculeBoundaryInfo.getZmin();
                    float tmpZmax = tmpMoleculeBoundaryInfo.getZmax();
                    boolean tmpIsMinSmallerMax = tmpZmin < tmpZmax;
                    boolean tmpIsMinGreaterMax = tmpZmin > tmpZmax;
                    for (int i = tmpMoleculeBoundaryInfo.getFirstIndex(); i < tmpMoleculeBoundaryInfo.getExclusiveLastIndex(); i++) {
                        if (tmpIsMinSmallerMax) {
                            if (tmpR_z[i] < tmpZmin) tmpR_z[i] = tmpZmin;
                            if (tmpR_z[i] > tmpZmax) tmpR_z[i] = tmpZmax;
                        } else if (tmpIsMinGreaterMax) {
                            if (tmpR_z[i] < tmpZmin && tmpR_z[i] > tmpZmax) {
                                if (tmpZmin - tmpR_z[i] <= tmpR_z[i] - tmpZmax) tmpR_z[i] = tmpZmin;
                                if (tmpZmin - tmpR_z[i] > tmpR_z[i] - tmpZmax) tmpR_z[i] = tmpZmax;
                            }
                        }
                        // Min = Max: Do nothing!
                    }
                }
            }
        }        
    }

    /**
     * Corrects molecule spheres if necessary.
     */
    private void correctMoleculeSpheres() {
        if (this.parameters.getChemicalSystemDescription().getMoleculeSphereInfos() != null) {
            float[] tmpR_x = this.parameters.getParticleArrays().getR_x();
            float[] tmpR_y = this.parameters.getParticleArrays().getR_y();
            float[] tmpR_z = this.parameters.getParticleArrays().getR_z();
            for (MoleculeSphereInfo tmpMoleculeSphereInfo: this.parameters.getChemicalSystemDescription().getMoleculeSphereInfos()) {
                int tmpFirstIndex = tmpMoleculeSphereInfo.getFirstIndex();
                int tmpExclusiveLastIndex = tmpMoleculeSphereInfo.getExclusiveLastIndex();
                for (int i = tmpFirstIndex; i < tmpExclusiveLastIndex; i++) {
                    float tmpDeltaX = tmpR_x[i] - tmpMoleculeSphereInfo.getSphereCenterX();
                    float tmpDeltaY = tmpR_y[i] - tmpMoleculeSphereInfo.getSphereCenterY();
                    float tmpDeltaZ = tmpR_z[i] - tmpMoleculeSphereInfo.getSphereCenterZ();
                    float tmpDistanceSquare = tmpDeltaX * tmpDeltaX + tmpDeltaY * tmpDeltaY + tmpDeltaZ * tmpDeltaZ;
                    boolean tmpIsInsideSphere = tmpDistanceSquare <= tmpMoleculeSphereInfo.getSphereRadiusSquare();
                    if (tmpMoleculeSphereInfo.isExclusiveSphere() && tmpIsInsideSphere) {
                        // <editor-fold defaultstate="collapsed" desc="Exclusive sphere (no particle allowed INSIDE sphere)">
                        float tmpDistance = (float) FastMath.sqrt(tmpDistanceSquare);
                        // NOTE: Unit vector tmpR_new_0 points from the center of the sphere in radial direction towards aR[i]
                        float tmpR_new_0_x = tmpDeltaX/tmpDistance;
                        float tmpR_new_0_y = tmpDeltaY/tmpDistance;
                        float tmpR_new_0_z = tmpDeltaZ/tmpDistance;
                        // New point outside sphere in radial direction
                        float tmpNewLength = tmpMoleculeSphereInfo.getSphereDiameter() - tmpDistance;
                        tmpR_x[i] = tmpMoleculeSphereInfo.getSphereCenterX() + tmpR_new_0_x * tmpNewLength;
                        tmpR_y[i] = tmpMoleculeSphereInfo.getSphereCenterY() + tmpR_new_0_y * tmpNewLength;
                        tmpR_z[i] = tmpMoleculeSphereInfo.getSphereCenterZ() + tmpR_new_0_z * tmpNewLength;
                        // </editor-fold>
                    } else if (!tmpMoleculeSphereInfo.isExclusiveSphere() && !tmpIsInsideSphere) {
                        // <editor-fold defaultstate="collapsed" desc="Inclusive sphere (no particle allowed OUTSIDE sphere)">
                        float tmpDistance = (float) FastMath.sqrt(tmpDistanceSquare);
                        // NOTE: Unit vector tmpR_new_0 points from the center of the sphere in radial direction towards aR[i]
                        float tmpR_new_0_x = tmpDeltaX/tmpDistance;
                        float tmpR_new_0_y = tmpDeltaY/tmpDistance;
                        float tmpR_new_0_z = tmpDeltaZ/tmpDistance;
                        // New point inside sphere in radial direction
                        float tmpNewLength = tmpMoleculeSphereInfo.getSphereDiameter() - tmpDistance;
                        tmpR_x[i] = tmpMoleculeSphereInfo.getSphereCenterX() + tmpR_new_0_x * tmpNewLength;
                        tmpR_y[i] = tmpMoleculeSphereInfo.getSphereCenterY() + tmpR_new_0_y * tmpNewLength;
                        tmpR_z[i] = tmpMoleculeSphereInfo.getSphereCenterZ() + tmpR_new_0_z * tmpNewLength;
                        // </editor-fold>
                    }
                }
            }
        }        
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Memory release related methods">
    /**
     * Releases memory and finishes logging
     */
    private void releaseMemoryAndFinishLogging() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("DpdSimulationTask.releaseMemoryAndFinishLogging", tmpId);
        // </editor-fold>
        // Do NOT set this.progressMonitor to null since this object comes from 
        // the calling method!
        this.factory = null;
        this.simulationOutput = null;
        this.parameters = null;
        this.particlePositionPool = null;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("DpdSimulationTask.releaseMemoryAndFinishLogging", tmpId);
        // </editor-fold>
        this.simulationLogger.finish();
        this.simulationLogger = null;
    }
    // </editor-fold>
    // </editor-fold>

}
