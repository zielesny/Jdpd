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
import de.gnwi.jdpd.movement.MoleculeSphereInfo;
import de.gnwi.jdpd.parameters.ChemicalSystemDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.parameters.SimulationDescription;
import de.gnwi.jdpd.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpd.utilities.GravitationalAcceleration;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Self-consistent Modified Velocity-Verlet (SCMVV) time step calculator
 * Literature:
 * Besold, G., Vattulainen, I. T., Karttunen, M., and Polson, J. M., 
 * Towards better integrators for dissipative particle dynamics simulations, 
 * Physical Review E. 62(6), 2000, 7611-7614
 * 
 * @author Achim Zielesny
 */
public class ScmvvTimeStepCalculator implements ITimeStepCalculator {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Factory for new objects
     */
    private final Factory factory;

    /**
     * Number of self-consistent iterations for Self-consistent Modified Velocity-Verlet (SCMVV) integration
     */
    private final int selfConsistentIterationNumber;
    
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
     * Full force accumulator
     */
    private final ForceAccumulator fullForceAccumulator;
    
    /**
     * Dissipative force accumulator
     */
    private final ForceAccumulator dissipativeForceAccumulator;
    
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
     * Molecule sphere infos
     */
    private final MoleculeSphereInfo[] moleculeSphereInfos;
    
    /**
     * True: Molecule sphere infos are defined, false: Otherwise
     */
    private final boolean hasMoleculeSphereInfos;
    
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
     * @param aSelfConsistentIterationNumber SelfConsistentIterationNumber parameter for Groot-Warren Modified Velocity-Verlet (GWMVV) integration
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
    public ScmvvTimeStepCalculator(
        Factory aFactory,
        int aSelfConsistentIterationNumber,
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
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aFactory is null.");
        }
        if (aSelfConsistentIterationNumber < 1) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aSelfConsistentIterationNumber < 1.");
        }
        if (aSimulationOutput == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aSimulationOutput is null.");
        }
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aSimulationLogger is null.");
        }
        if (aParameters == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aParameters is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aParallelizationInfo is null.");
        }
        if (aRandomNumberSeed == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aRandomNumberSeed is null.");
        }
        if (aParticlePositionPool == null) {
            throw new IllegalArgumentException("ScmvvTimeStepCalculator.Constructor: aParticlePositionPool is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("ScmvvTimeStepCalculator.Constructor", tmpId);
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
        this.moleculeSphereInfos = this.chemicalSystemDescription.getMoleculeSphereInfos();
        this.hasMoleculeSphereInfos = this.moleculeSphereInfos != null;
        this.moleculeVelocityFixationInfos = this.chemicalSystemDescription.getMoleculeVelocityFixationInfos();
        this.hasMoleculeVelocityFixationInfos = this.moleculeVelocityFixationInfos != null;
        this.moleculeAccelerationInfos = this.chemicalSystemDescription.getMoleculeAccelerationInfos();
        this.hasMoleculeAccelerationInfos = this.moleculeAccelerationInfos != null;

        this.gravitationalAcceleration = this.parameters.getInteractionDescription().getGravitationalAcceleration();

        this.numberOfInitialVelocityScalingSteps = this.parameters.getSimulationDescription().getNumberOfInitialVelocityScalingSteps();
        
        AtomicInteger tmpDummyRandomNumberSeed = new AtomicInteger(0);
        
        // NOTE: tmpParticlePairDpdConservativeRandomForceCalculator DOES use random numbers thus provide this.randomNumberSeed
        IParticlePairForceCalculator tmpParticlePairDpdForceFullCalculator =
            this.factory.getParticlePairScmvvDpdForceFullCalculator(
                this.simulationLogger,
                this.chemicalSystemDescription.getBoxSize(), 
                this.simulationDescription.getPeriodicBoundaries(), 
                this.factory.getDpdCutOffLength(),
                aParallelizationInfo, 
                aRandomNumberSeed
            );
        tmpParticlePairDpdForceFullCalculator.setParticlePairDistanceParametersCacheActivity(this.isCacheUsage);

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
        
        this.fullForceAccumulator = 
            new ForceAccumulator(
                this.simulationLogger,
                tmpParticlePairDpdForceFullCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                tmpBondForceCalculator,
                tmpParticlePairElectrostaticsForceCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
            );
        
        // NOTE: tmpParticlePairDpdDissipativeForceCalculator does NOT use any random numbers thus provide dummy
        IParticlePairForceCalculator tmpParticlePairDpdDissipativeForceCalculator =
            this.factory.getParticlePairScmvvDpdForceDissipativeCalculator(
                this.simulationLogger,
                this.chemicalSystemDescription.getBoxSize(), 
                this.simulationDescription.getPeriodicBoundaries(), 
                this.factory.getDpdCutOffLength(),
                aParallelizationInfo, 
                tmpDummyRandomNumberSeed
            );

        if (this.isCacheUsage) {
            this.dissipativeForceAccumulator = 
                new ForceAccumulator(
                    this.simulationLogger,
                    tmpParticlePairDpdDissipativeForceCalculator,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_CACHE,
                    null,
                    null,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.UNDEFINED
                );
        } else {
            // IMPORTANT: Dissipative force re-use particle cell assignments
            // thus calculation mode is WITHOUT_PARTICLE_CELL_ASSIGNMENTS
            this.dissipativeForceAccumulator = 
                new ForceAccumulator(
                    this.simulationLogger,
                    tmpParticlePairDpdDissipativeForceCalculator,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.WITHOUT_PARTICLE_CELL_ASSIGNMENTS,
                    null,
                    null,
                    ParticlePairInteractionCalculator.CellBasedCalculationMode.UNDEFINED
                );
        }
        
        this.particlePositionPool = aParticlePositionPool;

        this.selfConsistentIterationNumber = aSelfConsistentIterationNumber;
        
        IHarmonicBondPropertyCalculator tmpBondPotentialCalculator = null;
        if (this.fullForceAccumulator.getHarmonicBondForceConservativeCalculator() != null) {
            tmpBondPotentialCalculator = this.factory.getHarmonicBondPotentialCalculator(this.fullForceAccumulator.getHarmonicBondForceConservativeCalculator());
        }
        IParticlePairInteractionCalculator tmpParticlePairElectrostaticsPotentialCalculator = null;
        if (this.fullForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator() != null) {
            tmpParticlePairElectrostaticsPotentialCalculator = this.factory.getParticlePairElectrostaticsPotentialCalculator(this.fullForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator());
        }
        this.potentialAccumulator = 
            new PotentialAccumulator(
                this.simulationLogger,
                this.factory.getParticlePairDpdPotentialCalculator(this.fullForceAccumulator.getParticlePairDpdForceCalculator()),
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                tmpBondPotentialCalculator,
                tmpParticlePairElectrostaticsPotentialCalculator,
                ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
            );        

        if (this.simulationLogger.isLogLevel(ILogger.PARTICLE)) {
            IParticlePairForceCalculator tmpParticlePairDpdForceConservativeCalculator =
                this.factory.getParticlePairDpdForceConservativeCalculator(
                    this.simulationLogger,
                    this.chemicalSystemDescription.getBoxSize(), 
                    this.simulationDescription.getPeriodicBoundaries(), 
                    this.factory.getDpdCutOffLength(),
                    aParallelizationInfo, 
                    tmpDummyRandomNumberSeed
                );
            // IMPORTANT: Set cache usage
            tmpParticlePairDpdForceConservativeCalculator.setParticlePairDistanceParametersCacheActivity(true);
            
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
                    tmpParticlePairDpdForceConservativeCalculator,
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
                    this.factory.getParticlePairDpdForceConservativeCalculator(this.fullForceAccumulator.getParticlePairDpdForceCalculator());
                if (this.simulationDescription.isInitialPotentialEnergyMinimizationWithAllForces()) {
                    this.simulationLogger.appendIntermediateResults("Initial potential energy minimization with all forces");
                    Utils.calculateInitialPotentialEnergyMinimizationSteps(
                        this.simulationLogger,
                        this.parameters,
                        this.factory,
                        new ForceAccumulator(
                            this.simulationLogger,
                            tmpParticlePairDpdForceConservativeMinStepCalculator,
                            ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                            this.fullForceAccumulator.getHarmonicBondForceConservativeCalculator(),
                            this.fullForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator(),
                            ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
                        ),
                        this.potentialAccumulator,
                        this.simulationOutput,
                        this.particlePositionPool,
                        this.maximumNumberOfPositionCorrectionTrials
                    );
                } else {
                    this.simulationLogger.appendIntermediateResults("Initial potential energy minimization with DPD forces only");
                    Utils.calculateInitialPotentialEnergyMinimizationSteps(
                        this.simulationLogger,
                        this.parameters,
                        this.factory,
                        new ForceAccumulator(
                            this.simulationLogger,
                            tmpParticlePairDpdForceConservativeMinStepCalculator,
                            ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS,
                            null,
                            null,
                            ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_PARTICLE_CELL_ASSIGNMENTS
                        ),
                        this.potentialAccumulator,
                        this.simulationOutput,
                        this.particlePositionPool,
                        this.maximumNumberOfPositionCorrectionTrials
                    );
                }
                this.simulationOutput.setMinimizedParticlePositions(Utils.getParticlePositions(this.parameters, this.particlePositionPool));
                tmpParticlePairDpdForceConservativeMinStepCalculator.shutdownExecutorService();
            }
        }
        
        this.initialiseIntegration();
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("ScmvvTimeStepCalculator.Constructor", tmpId);
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
        this.simulationLogger.appendMethodCallStart("ScmvvTimeStepCalculator.calculate", tmpIdCalculate);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Initial call">
        if (this.isInitialCall) {
            // NOTE: this.fullForceAccumulator assigns conservative and random forces to f and dissipative forces to ftwo
            this.fullForceAccumulator.accumulate_f_and_fTwo(
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
        this.calculate_v_with_f_and_fOld(
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.particleArrays.getF_x(),
            this.particleArrays.getF_y(),
            this.particleArrays.getF_z(),
            this.particleArrays.getFtwo_x(),
            this.particleArrays.getFtwo_y(),
            this.particleArrays.getFtwo_z(),
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
                this.particleArrays.getRold_x(),
                this.particleArrays.getRold_y(),
                this.particleArrays.getRold_z()
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
        // Molecule spheres (AFTER Utils.correct_r_and_v_forMoleculeBoundaries())
        if (this.hasMoleculeSphereInfos) {
            Utils.correct_r_and_v_forMoleculeSpheres(
                aCurrentTimeStep,
                this.moleculeSphereInfos,
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z()
            );
        }
        // NOTE: this.fullForceAccumulator assigns conservative and random forces to f and dissipative forces to ftwo
        this.fullForceAccumulator.accumulate_f_and_fTwo(
            this.particleArrays.getBondChunkArraysList(),
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.parameters
        );
        if (this.isCacheUsage) {
            this.dissipativeForceAccumulator.getParticlePairDpdForceCalculator().setParticlePairDistanceParametersCache(
                this.fullForceAccumulator.getParticlePairDpdForceCalculator()
            );
        } else {
            // IMPORTANT: Set particle assignments to cell boxes in dissipative force calculator
            this.dissipativeForceAccumulator.getParticlePairDpdForceCalculator().setParticleCellAssignments(
                this.fullForceAccumulator.getParticlePairDpdForceCalculator()
            );
        }
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
        this.calculate_vNew_with_v_and_f(
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.particleArrays.getVnew_x(),
            this.particleArrays.getVnew_y(),
            this.particleArrays.getVnew_z(),
            this.particleArrays.getF_x(),
            this.particleArrays.getF_y(),
            this.particleArrays.getF_z(),
            this.particleArrays.getDpdMasses(),
            this.simulationDescription.getTimeStepLengthHalf(),
            this.simulationDescription.isDpdUnitMass()
        );
        for (int i = 0; i < this.selfConsistentIterationNumber; i++) {
            this.calculate_v_with_vNew_and_fold(
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z(),
                this.particleArrays.getVnew_x(),
                this.particleArrays.getVnew_y(),
                this.particleArrays.getVnew_z(),
                this.particleArrays.getFtwo_x(),
                this.particleArrays.getFtwo_y(),
                this.particleArrays.getFtwo_z(),
                this.particleArrays.getDpdMasses(),
                this.simulationDescription.getTimeStepLengthHalf(),
                this.simulationDescription.isDpdUnitMass()
            );
            // <editor-fold defaultstate="collapsed" desc="SCMVV logging preparation">
            if (this.simulationLogger.isLogLevel(ILogger.SCMVV)) {
                if (aCurrentTimeStep == 1 || aCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    this.simulationLogger.setScmvvInformationAccumulation(true);
                    this.simulationLogger.setScmvvFdissParticleIndexPairCounter(new AtomicInteger());
                    this.simulationLogger.setScmvvFdissRij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissAbsRij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissVij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissAbsVij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissRijVij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissAbsRijVij_x_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissVdotR_Adder(new DoubleAdder());
                    this.simulationLogger.setScmvvFdissGammaFactor_Adder(new DoubleAdder());
                } else {
                    this.simulationLogger.setScmvvInformationAccumulation(false);
                }
            }
            // </editor-fold>
            // NOTE: this.dissipativeForceAccumulator assigns dissipative forces to ftwo
            this.dissipativeForceAccumulator.accumulate_f_and_fTwo(
                null,
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.parameters
            );
            // <editor-fold defaultstate="collapsed" desc="SCMVV logging">
            if (this.simulationLogger.isLogLevel(ILogger.SCMVV)) {
                if (aCurrentTimeStep == 1 || aCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissParticleIndexPairCounter = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissParticleIndexPairCounter().get())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissRij_x sum                = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissRij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissAbsRij_x sum             = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissAbsRij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissVij_x sum                = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissVij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissAbsVij_x sum             = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissAbsVij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissRijVij_x sum             = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissRijVij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissAbsRijVij_x sum          = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissAbsRijVij_x_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissVdotR sum                = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissVdotR_Adder().sum())
                    );
                    this.simulationLogger.appendScmvv(
                        "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                        String.valueOf(i) + 
                        ") ScmvvFdissGammaFactor sum          = " + 
                        String.valueOf(this.simulationLogger.getScmvvFdissGammaFactor_Adder().sum())
                    );
                }
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Particle force magnitude logging for dissipative force">
            if (this.simulationLogger.isLogLevel(ILogger.PARTICLE)) {
                if (aCurrentTimeStep == 1 || aCurrentTimeStep%this.parameters.getSimulationDescription().getTimeStepFrequencyForOutput() == 0) {
                    double[] tmpMinMeanMaxDpdForceDissipativeMagnitude = 
                        Utils.calculateMinMeanMax(
                            this.particleArrays.getFtwo_x(),
                            this.particleArrays.getFtwo_y(),
                            this.particleArrays.getFtwo_z()
                        );
                    if (tmpMinMeanMaxDpdForceDissipativeMagnitude != null) {
                        this.simulationLogger.appendParticleForceMagnitude(
                            "ScmvvTimeStepCalculator.calculate: (SC iteration " + 
                            String.valueOf(i) + 
                            ") DPD dissipat. F (min/mean/max) = " + 
                            String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[0]) +
                            " / " +
                            String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[1]) +
                            " / " +
                            String.valueOf(tmpMinMeanMaxDpdForceDissipativeMagnitude[2])
                        );
                    }
                }
            }
            // </editor-fold>
        }
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
            this.simulationLogger.appendVelocityScaleFactor("ScmvvTimeStepCalculator.calculate, velocity scale factor WITH scaling = " + String.valueOf(tmpVelocityScaleFactor));
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
            this.simulationLogger.appendVelocityScaleFactor("ScmvvTimeStepCalculator.calculate, velocity scale factor WITH scaling due to molecule acceleration = " + String.valueOf(tmpVelocityScaleFactor));
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
            this.simulationLogger.appendVelocityScaleFactor("ScmvvTimeStepCalculator.calculate, velocity scale factor WITHOUT scaling = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        }
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("ScmvvTimeStepCalculator.calculate", tmpIdCalculate);
        // </editor-fold>
    }
    
    /**
     * Executor services shutdown
     */
    @Override
    public void shutdownExecutorServices() {
        this.fullForceAccumulator.shutdownExecutorService();
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
     * (No checks are performed)
     */
    private void initialiseIntegration() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("ScmvvTimeStepCalculator.initialiseIntegration", tmpId);
        // </editor-fold>
        if (this.hasMoleculeFixationInfos) {
            // Copy r to rOld
            Utils.copyToOld(
                this.particleArrays.getR_x(),
                this.particleArrays.getR_y(),
                this.particleArrays.getR_z(),
                this.particleArrays.getRold_x(),
                this.particleArrays.getRold_y(),
                this.particleArrays.getRold_z()
            );
        }
        this.isInitialCall = true;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("ScmvvTimeStepCalculator.initialiseIntegration", tmpId);
        // </editor-fold>
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Calculation methods">
    /**
     * Calculates v with f and ftwo
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aF_x Current x-components of particle forces
     * @param aF_y Current y-components of particle forces
     * @param aF_z Current z-components of particle forces
     * @param aFtwo_x Force two x-components of particles
     * @param aFtwo_y Force two y-components of particles
     * @param aFtwo_z Force two z-components of particles
     * @param aDpdMass DPD masses of particles
     * @param aTimeStepLengthHalf Half of time step length
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    private void calculate_v_with_f_and_fOld(
        double[] aV_x,
        double[] aV_y,
        double[] aV_z,
        double[] aF_x,
        double[] aF_y,
        double[] aF_z,
        double[] aFtwo_x,
        double[] aFtwo_y,
        double[] aFtwo_z,
        double[] aDpdMass,
        double aTimeStepLengthHalf,
        boolean anIsDpdUnitMass) {
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] += aTimeStepLengthHalf * (aF_x[i] + aFtwo_x[i]);
                aV_y[i] += aTimeStepLengthHalf * (aF_y[i] + aFtwo_y[i]);
                aV_z[i] += aTimeStepLengthHalf * (aF_z[i] + aFtwo_z[i]);
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] += aTimeStepLengthHalf * (aF_x[i] + aFtwo_x[i]) / aDpdMass[i];
                aV_y[i] += aTimeStepLengthHalf * (aF_y[i] + aFtwo_y[i]) / aDpdMass[i];
                aV_z[i] += aTimeStepLengthHalf * (aF_z[i] + aFtwo_z[i]) / aDpdMass[i];
            }
        }
    }

    /**
     * Calculates vNew with v and f
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aVnew_x New x-components of particle velocities (may be changed)
     * @param aVnew_y New y-components of particle velocities (may be changed)
     * @param aVnew_z New z-components of particle velocities (may be changed)
     * @param aF_x Current x-components of particle forces
     * @param aF_y Current y-components of particle forces
     * @param aF_z Current z-components of particle forces
     * @param aDpdMass DPD masses of particles
     * @param aTimeStepLengthHalf Half time step length
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    private void calculate_vNew_with_v_and_f(
        double[] aV_x,
        double[] aV_y,
        double[] aV_z,
        double[] aVnew_x,
        double[] aVnew_y,
        double[] aVnew_z,
        double[] aF_x,
        double[] aF_y,
        double[] aF_z,
        double[] aDpdMass,
        double aTimeStepLengthHalf,
        boolean anIsDpdUnitMass) {
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                aVnew_x[i] = aV_x[i] + aTimeStepLengthHalf * aF_x[i];
                aVnew_y[i] = aV_y[i] + aTimeStepLengthHalf * aF_y[i];
                aVnew_z[i] = aV_z[i] + aTimeStepLengthHalf * aF_z[i];
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                aVnew_x[i] = aV_x[i] + aTimeStepLengthHalf * aF_x[i] / aDpdMass[i];
                aVnew_y[i] = aV_y[i] + aTimeStepLengthHalf * aF_y[i] / aDpdMass[i];
                aVnew_z[i] = aV_z[i] + aTimeStepLengthHalf * aF_z[i] / aDpdMass[i];
            }
        }
    }

    /**
     * Calculates v with vNew and f
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aVnew_x New x-components of particle velocities (may be changed)
     * @param aVnew_y New y-components of particle velocities (may be changed)
     * @param aVnew_z New z-components of particle velocities (may be changed)
     * @param aFtwo_x Force two x-components of particles
     * @param aFtwo_y Force two y-components of particles
     * @param aFtwo_z Force two z-components of particles
     * @param aDpdMass DPD masses of particles
     * @param aTimeStepLengthHalf Half time step length
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    private void calculate_v_with_vNew_and_fold(
        double[] aV_x,
        double[] aV_y,
        double[] aV_z,
        double[] aVnew_x,
        double[] aVnew_y,
        double[] aVnew_z,
        double[] aFtwo_x,
        double[] aFtwo_y,
        double[] aFtwo_z,
        double[] aDpdMass,
        double aTimeStepLengthHalf,
        boolean anIsDpdUnitMass) {
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] = aVnew_x[i] + aTimeStepLengthHalf * aFtwo_x[i];
                aV_y[i] = aVnew_y[i] + aTimeStepLengthHalf * aFtwo_y[i];
                aV_z[i] = aVnew_z[i] + aTimeStepLengthHalf * aFtwo_z[i];
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] = aVnew_x[i] + aTimeStepLengthHalf * aFtwo_x[i] / aDpdMass[i];
                aV_y[i] = aVnew_y[i] + aTimeStepLengthHalf * aFtwo_y[i] / aDpdMass[i];
                aV_z[i] = aVnew_z[i] + aTimeStepLengthHalf * aFtwo_z[i] / aDpdMass[i];
            }
        }
    }
    // </editor-fold>
    // </editor-fold>

}
