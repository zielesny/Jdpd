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
package de.gnwi.jdpdsp.samples.integrationType;

import de.gnwi.jdpd.interfaces.ILogLevel;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.accumulators.ForceAccumulator;
import de.gnwi.jdpdsp.accumulators.ParticleForceMagnitudeAccumulator;
import de.gnwi.jdpdsp.accumulators.PotentialAccumulator;
import de.gnwi.jdpdsp.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IOutput;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.particlePosition.ParticlePositionPool;
import de.gnwi.jdpdsp.utilities.Factory;
import de.gnwi.jdpdsp.utilities.Utils;
import de.gnwi.jdpdsp.interfaces.ITimeStepCalculator;
import de.gnwi.jdpdsp.movement.MoleculeFixationInfo;
import de.gnwi.jdpdsp.movement.MoleculeVelocityFixationInfo;
import de.gnwi.jdpdsp.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpdsp.interfaces.IParticlePairInteractionPnhlnCalculator;
import de.gnwi.jdpdsp.interfaces.IRandom;
import de.gnwi.jdpdsp.movement.MoleculeAccelerationInfo;
import de.gnwi.jdpdsp.movement.MoleculeBoundaryInfo;
import de.gnwi.jdpdsp.movement.MoleculeSphereInfo;
import de.gnwi.jdpdsp.parameters.ChemicalSystemDescription;
import de.gnwi.jdpdsp.parameters.InteractionDescription;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.parameters.ParticleArrays;
import de.gnwi.jdpdsp.parameters.SimulationDescription;
import de.gnwi.jdpdsp.samples.interactions.ParticlePairInteractionCalculator;
import de.gnwi.jdpdsp.utilities.GravitationalAcceleration;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.math3.util.FastMath;

/**
 * Nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) time step 
 * calculator
 * Literature:
 * B. Leimkuhler, X. Shang, 
 * On the numerical treatment of dissipative particle dynamics and related 
 * systems, 
 * Journal of Computational Physics 280, 2015, 72–95
 * 
 * @author Achim Zielesny
 */
public class PnhlnTimeStepCalculator implements ITimeStepCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Numeric constants
     */
    private static final float ONE = 1.0f;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Factory for new objects
     */
    private final Factory factory;

    /**
     * Mu parameter for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
     */
    private final float mu;
    
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
     * Interaction description
     */
    private final InteractionDescription interactionDescription;

    /**
     * S1MVV particle pair velocity update calculator
     */
    private final IParticlePairInteractionPnhlnCalculator particlePairPnhlnVelocityUpdatePlusGCalculator;
    
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
     * Random number generator
     */
    private final IRandom randomNumberGenerator;
    
    /**
     * Number of initial velocity scaling steps
     */
    private final int numberOfInitialVelocityScalingSteps;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
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
     * @param aMu Mu parameter for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
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
    public PnhlnTimeStepCalculator(
        Factory aFactory,
        float aMu,
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
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aFactory is null.");
        }
        if (aSimulationOutput == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aSimulationOutput is null.");
        }
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aSimulationLogger is null.");
        }
        if (aParameters == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aParameters is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aParallelizationInfo is null.");
        }
        if (aRandomNumberSeed == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aRandomNumberSeed is null.");
        }
        if (aParticlePositionPool == null) {
            throw new IllegalArgumentException("PnhlnTimeStepCalculator.Constructor: aParticlePositionPool is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("PnhlnTimeStepCalculator.Constructor", tmpId);
        // </editor-fold>
        
        this.maximumNumberOfPositionCorrectionTrials = aMaximumNumberOfPositionCorrectionTrials;
        this.factory = aFactory;
        this.mu = aMu;
        this.simulationOutput = aSimulationOutput;
        this.isCacheUsage = anIsCacheUsage;
        
        this.parameters = aParameters;
        this.particleArrays = this.parameters.getParticleArrays();
        this.simulationDescription = this.parameters.getSimulationDescription();
        this.chemicalSystemDescription = this.parameters.getChemicalSystemDescription();
        this.interactionDescription = this.parameters.getInteractionDescription();
        // IMPORTANT: Set Ksi with DPD gamma
        this.interactionDescription.setPnhln_Ksi(this.interactionDescription.getDpdGamma());

        this.randomNumberGenerator = aFactory.getNewOrJumpedRandomNumberGenerator(aRandomNumberSeed.incrementAndGet());
        
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
            // NOTE: tmpParticlePairElectrostaticsForceCalculator does NOT use any random numbers thus provide dummy
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
        
        // NOTE: this.particlePairPnhlnVelocityUpdatePlusGCalculator does NOT use any random numbers thus provide dummy
        this.particlePairPnhlnVelocityUpdatePlusGCalculator =
            this.factory.getParticlePairPnhlnVelocityUpdatePlusGCalculator(
                this.simulationLogger,
                this.chemicalSystemDescription.getBoxSize(), 
                this.simulationDescription.getPeriodicBoundaries(), 
                this.factory.getDpdCutOffLength(),
                aParallelizationInfo, 
                tmpDummyRandomNumberSeed
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

        if (this.simulationLogger.isLogLevel(ILogLevel.PARTICLE)) {
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
                            this.conservativeForceAccumulator.getHarmonicBondForceConservativeCalculator(),
                            this.conservativeForceAccumulator.getParticlePairElectrostaticsForceConservativeCalculator(),
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
        this.simulationLogger.appendMethodCallEnd("PnhlnTimeStepCalculator.Constructor", tmpId);
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
        this.simulationLogger.appendMethodCallStart("PnhlnTimeStepCalculator.calculate", tmpIdCalculate);
        // </editor-fold>
        boolean tmpIsVelocityScalingForMoleculeAcceleration = false;
        // NOTE: Use this.simulationDescription.getTimeStepLengthHalf() to account
        // for factor 0.5
        Utils.calculate_r_with_v(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.simulationDescription.getTimeStepLengthHalf()
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
        // Update particle pair velocities plus G
        this.interactionDescription.setPnhln_G_Calculation(true);
        ParticlePairInteractionCalculator.CellBasedCalculationMode tmpCellBasedCalculationMode;
        if (this.isCacheUsage) {
            this.particlePairPnhlnVelocityUpdatePlusGCalculator.setParticlePairDistanceParametersCache(
                this.conservativeForceAccumulator.getParticlePairDpdForceCalculator()
            );
            tmpCellBasedCalculationMode = ParticlePairInteractionCalculator.CellBasedCalculationMode.WITH_CACHE;
        } else {
            // IMPORTANT: Set particle assignments to cell boxes in velocity update calculator
            this.particlePairPnhlnVelocityUpdatePlusGCalculator.setParticleCellAssignments(
                this.conservativeForceAccumulator.getParticlePairDpdForceCalculator()
            );
            tmpCellBasedCalculationMode = ParticlePairInteractionCalculator.CellBasedCalculationMode.WITHOUT_PARTICLE_CELL_ASSIGNMENTS;
        }
        // Update velocities
        this.particlePairPnhlnVelocityUpdatePlusGCalculator.calculateParticlePairInteractions(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.parameters,
            tmpCellBasedCalculationMode
        );
        // Update Ksi
        float tmpFactor1 = (float) FastMath.exp(-this.interactionDescription.getDpdGamma() * this.simulationDescription.getTimeStepLength());
        float tmpFactor2 = this.particlePairPnhlnVelocityUpdatePlusGCalculator.getAccumulatedPnhln_G_AddersSum() / this.mu * this.simulationDescription.getTimeStepLengthHalf();
        this.interactionDescription.setPnhln_Ksi(
            tmpFactor1 * (this.interactionDescription.getPnhln_Ksi() + tmpFactor2) + 
            this.randomNumberGenerator.nextGaussian() * (float) FastMath.sqrt(this.interactionDescription.getTemperature() / this.mu * (ONE - tmpFactor1 * tmpFactor1)) + 
            tmpFactor2
        );
        // Update particle pair velocities WITHOUT G
        this.interactionDescription.setPnhln_G_Calculation(false);
        // Update velocities
        this.particlePairPnhlnVelocityUpdatePlusGCalculator.calculateParticlePairInteractions(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.parameters,
            tmpCellBasedCalculationMode
        );
        // NOTE: Use this.simulationDescription.getTimeStepLengthHalf() to account
        // for factor 0.5
        Utils.calculate_r_with_v(
            this.particleArrays.getR_x(),
            this.particleArrays.getR_y(),
            this.particleArrays.getR_z(),
            this.particleArrays.getV_x(),
            this.particleArrays.getV_y(),
            this.particleArrays.getV_z(),
            this.simulationDescription.getTimeStepLengthHalf()
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
        final float tmpVelocityScaleFactor;
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
            this.simulationLogger.appendVelocityScaleFactor("PnhlnTimeStepCalculator.calculate, velocity scale factor WITH scaling = " + String.valueOf(tmpVelocityScaleFactor));
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
            this.simulationLogger.appendVelocityScaleFactor("PnhlnTimeStepCalculator.calculate, velocity scale factor WITH scaling due to molecule acceleration = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        } else if (this.simulationLogger.isLogLevel(ILogLevel.V_SCALE)) {
            tmpVelocityScaleFactor = Utils.getVelocityScaleFactor(
                this.particleArrays.getV_x(),
                this.particleArrays.getV_y(),
                this.particleArrays.getV_z(),
                this.particleArrays.getDpdMasses(),
                this.parameters.getInteractionDescription().getTemperature(),
                this.simulationDescription.isDpdUnitMass()
            );
            // <editor-fold defaultstate="collapsed" desc="Velocity scale factor logging">
            this.simulationLogger.appendVelocityScaleFactor("PnhlnTimeStepCalculator.calculate, velocity scale factor WITHOUT scaling = " + String.valueOf(tmpVelocityScaleFactor));
            // </editor-fold>
        }
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("PnhlnTimeStepCalculator.calculate", tmpIdCalculate);
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
     * (No checks are performed)
     */
    private void initialiseIntegration() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("PnhlnTimeStepCalculator.initialiseIntegration", tmpId);
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
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("PnhlnTimeStepCalculator.initialiseIntegration", tmpId);
        // </editor-fold>
    }
    // </editor-fold>
    // </editor-fold>

}
