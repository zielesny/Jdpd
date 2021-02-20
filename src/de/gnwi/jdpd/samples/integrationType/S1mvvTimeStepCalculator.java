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
package de.gnwi.jdpd.samples.integrationType;

import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.accumulators.ForceAccumulator;
import de.gnwi.jdpd.accumulators.ParticleForceMagnitudeAccumulator;
import de.gnwi.jdpd.accumulators.PotentialAccumulator;
import de.gnwi.jdpd.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IOutput;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.particlePosition.ParticlePositionPool;
import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.utilities.Utils;
import de.gnwi.jdpd.interfaces.ITimeStepCalculator;
import de.gnwi.jdpd.movement.MoleculeFixationInfo;
import de.gnwi.jdpd.movement.MoleculeVelocityFixationInfo;
import de.gnwi.jdpd.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.movement.MoleculeAccelerationInfo;
import de.gnwi.jdpd.movement.MoleculeBoundaryInfo;
import de.gnwi.jdpd.parameters.ChemicalSystemDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.parameters.SimulationDescription;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.GravitationalAcceleration;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Shardlow S1 Modified Velocity-Verlet (S1MVV) time step calculator
 * Literature:
 * T. Shardlow, Splitting for Dissipative Particle Dynamics, 
 * SIAM J. Sci. Comp. 24 (4), 2003, 1267-1282
 * P. Nikunen, M. Karttunen, I. Vattulainen, 
 * How would you integrate the equations of motion in dissipative particle 
 * dynamics simulations?, 
 * Computer Physics Communications 153, 2003, 407–423
 * 
 * @author Achim Zielesny
 */
public class S1mvvTimeStepCalculator implements ITimeStepCalculator {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Factory for new objects
     */
    private final Factory factory;
    
    /**
     * Flag for cache usage, True: Cache will be used, false: Otherwise
     */
    private final boolean isCacheUsage;

    /**
     * Simulation output
     */
    private final IOutput simulationOutput;
    
    /**
     * Simulation logger
     */
    private final ILogger simulationLogger;
    
    /** 
     * Simulation parameters: Parameters
     */
    private final Parameters parameters;
    
    /**
     * Particle arrays
     */
    private final ParticleArrays particleArrays;
    
    /**
     * Simulation description
     */
    private final SimulationDescription simulationDescription;

    /**
     * Chemical system description
     */
    private final ChemicalSystemDescription chemicalSystemDescription;

    /**
     * S1MVV particle pair velocity update calculator
     */
    private final IParticlePairInteractionCalculator particlePairS1mvvVelocityUpdateCalculator;
    
    /**
     * Conservative force accumulator
     */
    private final ForceAccumulator conservativeForceAccumulator;
    
    /**
     * Potential accumulator
     */
    private final PotentialAccumulator potentialAccumulator;

    /**
     * Particle force magnitude accumulator
     */
    private final ParticleForceMagnitudeAccumulator particleForceMagnitudeAccumulator;
    
    /**
     * Particle position pools
     */
    private final ParticlePositionPool particlePositionPool;
    
    /**
     * Molecule fixation infos
     */
    private final MoleculeFixationInfo[] moleculeFixationInfos;
    
    /**
     * True: Molecule fixation infos are defined, false: Otherwise
     */
    private final boolean hasMoleculeFixationInfos;
    
    /**
     * Molecule boundary infos
     */
    private final MoleculeBoundaryInfo[] moleculeBoundaryInfos;
    
    /**
     * True: Molecule boundary infos are defined, false: Otherwise
     */
    private final boolean hasMoleculeBoundaryInfos;
    
    /**
     * Molecule velocity fixation infos
     */
    private final MoleculeVelocityFixationInfo[] moleculeVelocityFixationInfos;
    
    /**
     * True: Molecule velocity fixation infos are defined, false: Otherwise
     */
    private final boolean hasMoleculeVelocityFixationInfos;
    
    /**
     * Molecule acceleration infos
     */
    private final MoleculeAccelerationInfo[] moleculeAccelerationInfos;
    
    /**
     * True: Molecule acceleration infos are defined, false: Otherwise
     */
    private final boolean hasMoleculeAccelerationInfos;
    
    /**
     * Gravitational acceleration
     */
    private final GravitationalAcceleration gravitationalAcceleration;
    
    /**
     * Number of initial velocity scaling steps
     */
    private final int numberOfInitialVelocityScalingSteps;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * True: Initial call, false: Otherwise
     */
    private boolean isInitialCall;

    /**
     * Maximum number of position correction trials
     */
    private int maximumNumberOfPositionCorrectionTrials;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aFactory Factory for new objects
     * @param anIsCacheUsage Flag for cache usage, True: Cache will be used, false: Otherwise
     * @param aSimulationOutput Simulation output
     * @param aSimulationLogger Simulation logger
     * @param aParameters Parameters instance
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @param aParticlePositionPool Particle position pool (NOT allowed to be 
     * null)
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     * @throws IllegalArgumentException Thrown if argument is invalid
     */
    public S1mvvTimeStepCalculator(
        Factory aFactory,
        boolean anIsCacheUsage,
        IOutput aSimulationOutput,
        ILogger aSimulationLogger, 
        Parameters aParameters, 
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed,
        ParticlePositionPool aParticlePositionPool,
        int aMaximumNumberOfPositionCorrectionTrials
    ) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aFactory == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aFactory is null.");
        }
        if (aSimulationOutput == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aSimulationOutput is null.");
        }
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aSimulationLogger is null.");
        }
        if (aParameters == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aParameters is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aParallelizationInfo is null.");
        }
        if (aRandomNumberSeed == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aRandomNumberSeed is null.");
        }
        if (aParticlePositionPool == null) {
            throw new IllegalArgumentException("S1mvvTimeStepCalculator.Constructor: aParticlePositionPool is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("S1mvvTimeStepCalculator.Constructor", tmpId);
        // </editor-fold>
        
        this.maximumNumberOfPositionCorrectionTrials = aMaximumNumberOfPositionCorrectionTrials;
        this.factory = aFactory;
        this.simulationOutput = aSimulationOutput;
        this.isCacheUsage = anIsCacheUsage;
        
        this.parameters = aParameters;
        this.particleArrays = this.parameters.getParticleArrays();
        this.simulationDescription = this.parameters.getSimulationDescription();
        this.chemicalSystemDescription = this.parameters.getChemicalSystemDescription();
        
        this.moleculeFixationInfos = this.chemicalSystemDescription.getMoleculeFixationInfos();
        this.hasMoleculeFixationInfos = this.moleculeFixationInfos != null;
        this.moleculeBoundaryInfos = this.chemicalSystemDescription.getMoleculeBoundaryInfos();
        this.hasMoleculeBoundaryInfos = this.moleculeBoundaryInfos != null;
        this.moleculeVelocityFixationInfos = this.chemicalSystemDescription.getMoleculeVelocityFixationInfos();
        this.hasMoleculeVelocityFixationInfos = this.moleculeVelocityFixationInfos != null;
        this.moleculeAccelerationInfos = this.chemicalSystemDescription.getMoleculeAccelerationInfos();
        this.hasMoleculeAccelerationInfos = this.moleculeAccelerationInfos != null;

        this.gravitationalAcceleration = this.parameters.getInteractionDescription().getGravitationalAcceleration();

        this.numberOfInitialVelocityScalingSteps = this.parameters.getSimulationDescription().getNumberOfInitialVelocityScalingSteps();

        AtomicInteger tmpDummyRandomNumberSeed = new AtomicInteger(0);
        
        // NOTE: tmpParticlePairDpdForceConservativeCalculator does NOT use any random numbers thus provide dummy
        IParticlePairForceCalculator tmpParticlePairDpdForceConservativeCalculator =
            this.factory.getParticlePairDpdForceConservativeCalculator(
                this.simulationLogger,
                this.chemicalSystemDescription.getBoxSize(), 
                this.simulationDescription.getPeriodicBoundaries(), 
                this.factory.getDpdCutOffLength(),
                aParallelizationInfo, 
                tmpDummyRandomNumberSeed
            );
        tmpParticlePairDpdForceConservativeCalculator.setParticlePairDistanceParametersCacheActivity(this.isCacheUsage);
        
        IHarmonicBondForceCalculator tmpBondForceCalculator = null;
        if (this.particleArrays.hasBonds()) {
            tmpBondForceCalculator = 
                this.factory.getHarmonicBondForceConservativeCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    aParallelizationInfo
                );
        }
        
        IParticlePairForceCalculator tmpParticlePairElectrostaticsForceCalculator = null;
        if (this.particleArrays.hasChargedParticles() && this.parameters.getInteractionDescription().hasElectrostatics()) {
            // NOTE: ParticlePairAdHocElectrostaticsForceCalculator does NOT use any random numbers thus provide dummy
            tmpParticlePairElectrostaticsForceCalculator =
                this.factory.getParticlePairElectrostaticsForceConservativeCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    this.parameters.getInteractionDescription().getElectrostatics().getCutOffLength(),
                    aParallelizationInfo, 
                    tmpDummyRandomNumberSeed
                );
        }
        
        this.conservativeForceAccumulator = 
            new ForceAccumulator(
                this.simulationLogger,
                tmpParticlePairDpdForceConservativeCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                tmpBondForceCalculator,
                tmpParticlePairElectrostaticsForceCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
            );
        
        // NOTE: this.particlePairS1mvvVelocityUpdateCalculator DOES use random numbers thus provide aRandomNumberSeed
        this.particlePairS1mvvVelocityUpdateCalculator =
            this.factory.getParticlePairS1mvvVelocityUpdateCalculator(
                this.simulationLogger,
                this.chemicalSystemDescription.getBoxSize(), 
                this.simulationDescription.getPeriodicBoundaries(), 
                this.factory.getDpdCutOffLength(),
                aParallelizationInfo, 
                aRandomNumberSeed
            );

        this.particlePositionPool = aParticlePositionPool;

        IHarmonicBondPropertyCalculator tmpBondPotentialCalculator = null;
        if (this.conservativeForceAccumulator.getHarmonicBondForceConservativeCalculator() != null) {
            tmpBondPotentialCalculator = this.factory.getHarmonicBondPotentialCalculator(this.conservativeForceAccumulator.getHarmonicBondForceConservativeCalculator());
        }
        IParticlePairInteractionCalculator tmpParticlePairElectrostaticsPotentialCalculator = null;
        if (this.conservativeForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator() != null) {
            tmpParticlePairElectrostaticsPotentialCalculator = this.factory.getParticlePairElectrostaticsPotentialCalculator(this.conservativeForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator());
        }
        this.potentialAccumulator = 
            new PotentialAccumulator(
                this.simulationLogger,
                this.factory.getParticlePairDpdPotentialCalculator(this.conservativeForceAccumulator.getParticlePairDpdForceCalculator()),
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                tmpBondPotentialCalculator,
                tmpParticlePairElectrostaticsPotentialCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
            );        

        if (this.simulationLogger.isLogLevel(ILogger.PARTICLE)) {
            IParticlePairForceCalculator tmpParticlePairDpdForceConservativeCalculator2 =
                this.factory.getParticlePairDpdForceConservativeCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    this.factory.getDpdCutOffLength(),
                    aParallelizationInfo, 
                    tmpDummyRandomNumberSeed
                );
            // IMPORTANT: Set cache usage
            tmpParticlePairDpdForceConservativeCalculator2.setParticlePairDistanceParametersCacheActivity(true);
            
            // NOTE: tmpParticlePairDpdFullForceCalculator DOES use random numbers thus provide random number seed
            IParticlePairForceCalculator tmpParticlePairDpdForceRandomCalculator =
                this.factory.getParticlePairDpdForceRandomCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    this.factory.getDpdCutOffLength(),
                    aParallelizationInfo, 
                    aRandomNumberSeed
                );
            IParticlePairForceCalculator tmpParticlePairDpdForceDissipativeCalculator =
                this.factory.getParticlePairDpdForceDissipativeCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    this.factory.getDpdCutOffLength(),
                    aParallelizationInfo, 
                    tmpDummyRandomNumberSeed
                );

            this.particleForceMagnitudeAccumulator = 
                new ParticleForceMagnitudeAccumulator(
                    this.simulationLogger,
                    tmpParticlePairDpdForceConservativeCalculator2,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                    tmpParticlePairDpdForceRandomCalculator,
                    tmpParticlePairDpdForceDissipativeCalculator,    
                    tmpBondForceCalculator,
                    tmpParticlePairElectrostaticsForceCalculator,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
                );
        } else {
            this.particleForceMagnitudeAccumulator = null;
        }
        
        if (!this.parameters.hasRestartInfo()) {
            this.simulationOutput.setStartParticlePositions(Utils.getParticlePositions(this.parameters, this.particlePositionPool));
            if (this.simulationDescription.getInitialPotentialEnergyMinimizationStepNumber() > 0) {
                IParticlePairForceCalculator tmpParticlePairDpdForceConservativeMinStepCalculator = 
                    this.factory.getParticlePairDpdForceConservativeCalculator(this.conservativeForceAccumulator.getParticlePairDpdForceCalculator());
                Utils.calculateInitialPotentialEnergyMinimizationSteps(
                    this.simulationLogger,
                    this.parameters,
                    this.factory,
                    new ForceAccumulator(
                        this.simulationLogger,
                        tmpParticlePairDpdForceConservativeMinStepCalculator,
                        ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                        this.conservativeForceAccumulator.getHarmonicBondForceConservativeCalculator(),
                        this.conservativeForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator(),
                        ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
                    ),
                    this.potentialAccumulator,
                    this.simulationOutput,
                    this.particlePositionPool,
                    this.maximumNumberOfPositionCorrectionTrials
                );
                this.simulationOutput.setMinimizedParticlePositions(Utils.getParticlePositions(this.parameters, this.particlePositionPool));
                tmpParticlePairDpdForceConservativeMinStepCalculator.shutdownExecutorService();
            }
        }
        
        this.initialiseIntegration();
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("S1mvvTimeStepCalculator.Constructor", tmpId);
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Calculates single time step
     * 
     * @param aCurrentTimeStep Current time step
     */
    @Override
    public void calculate(int aCurrentTimeStep) {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpIdCalculate = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("S1mvvTimeStepCalculator.calculate", tmpIdCalculate);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Initial call">
        if (this.isInitialCall) {
            this.conservativeForceAccumulator.accumulate_f_and_fTwo(
                this.particleArrays.getBondChunkArraysList(),
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.parameters
            );
            if (this.gravitationalAcceleration.isGravitationalAcceleration()) {
                Utils.addGravitationalForceTo_f(
                    this.gravitationalAcceleration, 
                    this.particleArrays.getF_x(),
                    this.particleArrays.getF_y(),
                    this.particleArrays.getF_z(),
                    this.particleArrays.getDpdMasses(),
                    this.simulationDescription.isDpdUnitMass()
                );
            }
            if (this.hasMoleculeAccelerationInfos) {
                for (MoleculeAccelerationInfo tmpMoleculeAccelerationInfo : this.moleculeAccelerationInfos) {
                    if (aCurrentTimeStep > 1 && aCurrentTimeStep <= tmpMoleculeAccelerationInfo.getMaxTimeStep()  && aCurrentTimeStep%tmpMoleculeAccelerationInfo.getFrequency() == 0) {
                        Utils.addMoleculeAccelerationTo_f(
                            tmpMoleculeAccelerationInfo,
                            this.particleArrays.getF_x(),
                            this.particleArrays.getF_y(),
                            this.particleArrays.getF_z(),
                            this.particleArrays.getDpdMasses(),
                            this.simulationDescription.isDpdUnitMass()
                        );
                    }
                }
            }
            this.isInitialCall = false;
        }
        // </editor-fold>
        boolean tmpIsVelocityScalingForMoleculeAcceleration = false;
        ParticlePairInteractionCalculator.CellBasedCalculationMode tmpCellBasedCalculationMode;
        if (this.isCacheUsage) {
            this.particlePairS1mvvVelocityUpdateCalculator.setParticlePairDistanceParametersCache(
                this.conservativeForceAccumulator.getParticlePairDpdForceCalculator()
            );
            tmpCellBasedCalculationMode = ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_CACHE;
        } else {
            // IMPORTANT: Set particle assignments to cell boxes in velocity update calculator
            this.particlePairS1mvvVelocityUpdateCalculator.setParticleCellAssignments(
                this.conservativeForceAccumulator.getParticlePairDpdForceCalculator()
            );
            tmpCellBasedCalculationMode = ParticlePairInteractionCalculator.CellBasedCalculationMode.WITHOUT_PARTICLE_CELL_ASSIGNMENTS;
        }
        
        // Update velocities
        this.particlePairS1mvvVelocityUpdateCalculator.calculateParticlePairInteractions(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.parameters,
            tmpCellBasedCalculationMode
        );
        
        Utils.calculate_v_with_f(
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.particleArrays.getF_x(),
            this.particleArrays.getF_y(),
            this.particleArrays.getF_z(),
            this.particleArrays.getDpdMasses(),
            this.simulationDescription.getTimeStepLengthHalf(),
            this.simulationDescription.isDpdUnitMass()
        );
        Utils.calculate_r_with_v(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.simulationDescription.getTimeStepLength()
        );
        // Molecule fixations
        if (this.hasMoleculeFixationInfos) {
            Utils.copy_rOld_to_r_forFixedMolecules(
                aCurrentTimeStep,
                this.moleculeFixationInfos,
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.particleArrays.getROld_x(),
                this.particleArrays.getROld_y(),
                this.particleArrays.getROld_z()
            );
        }
        Utils.correct_r_and_v(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.chemicalSystemDescription.getBoxSize(),
            this.simulationDescription.getPeriodicBoundaries(),
            this.maximumNumberOfPositionCorrectionTrials
        );
        // Molecule boundaries (AFTER Utils.correct_r_and_v())
        if (this.hasMoleculeBoundaryInfos) {
            Utils.correct_r_and_v_forMoleculeBoundaries(
                aCurrentTimeStep,
                this.moleculeBoundaryInfos,
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z()
            );
        }
        this.conservativeForceAccumulator.accumulate_f_and_fTwo(
            this.particleArrays.getBondChunkArraysList(),
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.parameters
        );
        if (this.gravitationalAcceleration.isGravitationalAcceleration()) {
            Utils.addGravitationalForceTo_f(
                this.gravitationalAcceleration, 
                this.particleArrays.getF_x(),
                this.particleArrays.getF_y(),
                this.particleArrays.getF_z(),
                this.particleArrays.getDpdMasses(),
                this.simulationDescription.isDpdUnitMass()
            );
        }
        if (this.hasMoleculeAccelerationInfos) {
            for (MoleculeAccelerationInfo tmpMoleculeAccelerationInfo : this.moleculeAccelerationInfos) {
                if (aCurrentTimeStep > 1 && aCurrentTimeStep <= tmpMoleculeAccelerationInfo.getMaxTimeStep()  && aCurrentTimeStep%tmpMoleculeAccelerationInfo.getFrequency() == 0) {
                    Utils.addMoleculeAccelerationTo_f(
                        tmpMoleculeAccelerationInfo,
                        this.particleArrays.getF_x(),
                        this.particleArrays.getF_y(),
                        this.particleArrays.getF_z(),
                        this.particleArrays.getDpdMasses(),
                        this.simulationDescription.isDpdUnitMass()
                    );
                    tmpIsVelocityScalingForMoleculeAcceleration = true;
                }
            }
        }
        Utils.calculate_v_with_f(
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.particleArrays.getF_x(),
            this.particleArrays.getF_y(),
            this.particleArrays.getF_z(),
            this.particleArrays.getDpdMasses(),
            this.simulationDescription.getTimeStepLengthHalf(),
            this.simulationDescription.isDpdUnitMass()
        );
        // Molecule velocity fixations
        if (this.hasMoleculeVelocityFixationInfos) {
            Utils.fill_v_forFixedMolecules(
                aCurrentTimeStep,
                this.moleculeVelocityFixationInfos,
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z()
            );            
        }
        final double tmpVelocityScaleFactor;
        if (aCurrentTimeStep <= this.numberOfInitialVelocityScalingSteps) {
            tmpVelocityScaleFactor = Utils.scale_v(
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z(),
                this.particleArrays.getDpdMasses(),
                this.parameters.getInteractionDescription().getTemperature(),
                this.simulationDescription.isDpdUnitMass()
            );
            // <editor-fold defaultstate="collapsed" desc="Velocity scale factor logging">
            this.simulationLogger.appendVelocityScaleFactor("S1mvvTimeStepCalculator.calculate, velocity scale factor WITH scaling = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        } else if (tmpIsVelocityScalingForMoleculeAcceleration) {
            tmpVelocityScaleFactor = Utils.scale_v(
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z(),
                this.particleArrays.getDpdMasses(),
                this.parameters.getInteractionDescription().getTemperature(),
                this.simulationDescription.isDpdUnitMass()
            );
            // <editor-fold defaultstate="collapsed" desc="Velocity scale factor logging">
            this.simulationLogger.appendVelocityScaleFactor("S1mvvTimeStepCalculator.calculate, velocity scale factor WITH scaling due to molecule acceleration = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        } else if (this.simulationLogger.isLogLevel(ILogger.V_SCALE)) {
            tmpVelocityScaleFactor = Utils.getVelocityScaleFactor(
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z(),
                this.particleArrays.getDpdMasses(),
                this.parameters.getInteractionDescription().getTemperature(),
                this.simulationDescription.isDpdUnitMass()
            );
            // <editor-fold defaultstate="collapsed" desc="Velocity scale factor logging">
            this.simulationLogger.appendVelocityScaleFactor("S1mvvTimeStepCalculator.calculate, velocity scale factor WITHOUT scaling = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        }
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("S1mvvTimeStepCalculator.calculate", tmpIdCalculate);
        // </editor-fold>
    }
    
    /**
     * Executor services shutdown
     */
    @Override
    public void shutdownExecutorServices() {
        this.conservativeForceAccumulator.shutdownExecutorService();
        this.potentialAccumulator.shutdownExecutorService();
        if (this.particleForceMagnitudeAccumulator != null) {
            this.particleForceMagnitudeAccumulator.shutdownExecutorService();
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Potential accumulator
     * 
     * @return Potential accumulator
     */
    @Override
    public PotentialAccumulator getPotentialAccumulator() {
        return this.potentialAccumulator;
    }
    
    /**
     * Particle force magnitude accumulator
     * 
     * @return Particle force magnitude accumulator
     */
    @Override
    public ParticleForceMagnitudeAccumulator getParticleForceMagnitudeAccumulator() {
        return this.particleForceMagnitudeAccumulator;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    // <editor-fold defaultstate="collapsed" desc="- Initialisation method">
    /**
     * Initialises integration
     * NOTE: No checks are performed.
     */
    private void initialiseIntegration() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("S1mvvTimeStepCalculator.initialiseIntegration", tmpId);
        // </editor-fold>
        if (this.hasMoleculeFixationInfos) {
            // Copy r to rOld
            Utils.copyToOld(
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.particleArrays.getROld_x(),
                this.particleArrays.getROld_y(),
                this.particleArrays.getROld_z()
            );
        }
        this.isInitialCall = true;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("S1mvvTimeStepCalculator.initialiseIntegration", tmpId);
        // </editor-fold>
    }
    // </editor-fold>
    // </editor-fold>

}
