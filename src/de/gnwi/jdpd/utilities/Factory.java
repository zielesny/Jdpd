/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2018  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.utilities;


import de.gnwi.jdpd.interfaces.IHarmonicBondForceCalculator;
import de.gnwi.jdpd.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IOutput;
import de.gnwi.jdpd.interfaces.IParticlePairForceCalculator;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionPnhlnCalculator;
import de.gnwi.jdpd.interfaces.IRandom;
import de.gnwi.jdpd.interfaces.ITimeStepCalculator;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.particlePosition.ParticlePositionPool;
import de.gnwi.jdpd.samples.harmonicBonds.HarmonicBondForceConservativeCalculator;
import de.gnwi.jdpd.samples.integrationType.GwmvvTimeStepCalculator;
import de.gnwi.jdpd.samples.harmonicBonds.HarmonicBondPotentialCalculator;
import de.gnwi.jdpd.samples.integrationType.PnhlnTimeStepCalculator;
import de.gnwi.jdpd.samples.integrationType.S1mvvTimeStepCalculator;
import de.gnwi.jdpd.samples.integrationType.ScmvvTimeStepCalculator;
import de.gnwi.jdpd.samples.interactions.electrostatics.ParticlePairElectrostaticsAdHocForceConservativeCalculator;
import de.gnwi.jdpd.samples.interactions.electrostatics.ParticlePairElectrostaticsAdHocPotentialCalculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairDpdForceConservativeCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairScmvvDpdForceDissipativeCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairGwmvvDpdForceFullCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairScmvvDpdForceFullCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairDpdPotentialCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.dpdCutoff1.ParticlePairS1mvvVelocityUpdateCutoff1Calculator;
import de.gnwi.jdpd.samples.interactions.nearestNeighbor.ParticlePairNearestNeighborCalculator;
import de.gnwi.jdpd.samples.random.ApacheCommonsRandom;
import de.gnwi.jdpd.samples.random.PcgRandom;
import java.util.concurrent.atomic.AtomicInteger;
import org.apache.commons.rng.simple.RandomSource;

/**
 * Factory for specific object instance creation
 * 
 * @author Achim Zielesny
 */
public class Factory {
    
    // <editor-fold defaultstate="collapsed" desc="Public Enums">
    /**
     * Random number generator type
     */
    public enum RandomType {

        /**
         * Apache Commons RNG (ACRNG)
         */
        ACRNG_ISAAC,
        ACRNG_JDK,
        ACRNG_KISS,
        ACRNG_MT,
        ACRNG_MT_64,
        ACRNG_MWC_256,
        ACRNG_SPLIT_MIX_64,
        ACRNG_TWO_CMRES,
        ACRNG_WELL_1024_A,
        ACRNG_WELL_19937_A,
        ACRNG_WELL_19937_C,
        ACRNG_WELL_44497_A,
        ACRNG_WELL_44497_B,
        ACRNG_WELL_512_A,
        ACRNG_XOR_SHIFT_1024_S,
        /**
         * PCG RNG
         */
        PCG_32;

        /**
         * Random number generator type representations
         * 
         * @return Random number generator type representations
         */
        public static String[] getRandomNumberGeneratorTypeRepresentations() {
            return new String[] {
                RandomType.ACRNG_ISAAC.toString(),
                RandomType.ACRNG_JDK.toString(),
                RandomType.ACRNG_KISS.toString(),
                RandomType.ACRNG_MT.toString(),
                RandomType.ACRNG_MT_64.toString(),
                RandomType.ACRNG_MWC_256.toString(),
                RandomType.ACRNG_SPLIT_MIX_64.toString(),
                RandomType.ACRNG_TWO_CMRES.toString(),
                RandomType.ACRNG_WELL_1024_A.toString(),
                RandomType.ACRNG_WELL_19937_A.toString(),
                RandomType.ACRNG_WELL_19937_C.toString(),
                RandomType.ACRNG_WELL_44497_A.toString(),
                RandomType.ACRNG_WELL_44497_B.toString(),
                RandomType.ACRNG_WELL_512_A.toString(),
                RandomType.ACRNG_XOR_SHIFT_1024_S.toString(),
                RandomType.PCG_32.toString()
            };
        }

        /**
         * Default random number generator type representation
         * 
         * @return Default random number generator type representation
         */
        public static String getDefaultRandomNumberGeneratorTypeRepresentation() {
            return Factory.RandomType.getDefaultRandomNumberGeneratorType().toString();
        }

        /**
         * Default random number generator type
         * 
         * @return Default random number generator type
         */
        public static RandomType getDefaultRandomNumberGeneratorType() {
            return RandomType.ACRNG_MT;
        }

        /**
         * Checks if aRandomNumberGeneratorTypeRepresentation is defined
         * 
         * @param aRandomNumberGeneratorTypeRepresentation
         * @return True: aRandomNumberGeneratorTypeRepresentation is defined, 
         * false: Otherwise
         */
        public static boolean isDefinedRandomNumberGeneratorTypeRepresentation(String aRandomNumberGeneratorTypeRepresentation) {
            // <editor-fold defaultstate="collapsed" desc="Checks">
            if (aRandomNumberGeneratorTypeRepresentation == null || aRandomNumberGeneratorTypeRepresentation.isEmpty()) {
                return false;
            }
            // </editor-fold>
            String[] tmpDefinedRandomNumberGeneratorTypes = Factory.RandomType.getRandomNumberGeneratorTypeRepresentations();
            for (String tmpDefinedRandomNumberGeneratorType : tmpDefinedRandomNumberGeneratorTypes) {
                if (aRandomNumberGeneratorTypeRepresentation.equals(tmpDefinedRandomNumberGeneratorType)) {
                    return true;
                }
            }
            return false;
        }

        /**
         * Returns RandomType of aRandomNumberGeneratorTypeRepresentation
         * 
         * @param aRandomNumberGeneratorTypeRepresentation
         * @return RandomType of aRandomNumberGeneratorTypeRepresentation or null if
         * not available
         */
        public static RandomType getRandomNumberGeneratorType(String aRandomNumberGeneratorTypeRepresentation) {
            Factory.RandomType tmpRandomType = null;
            if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_ISAAC.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_ISAAC;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_JDK.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_JDK;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_KISS.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_KISS;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_MT.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_MT;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_MT_64.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_MT_64;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_MWC_256.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_MWC_256;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_SPLIT_MIX_64.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_SPLIT_MIX_64;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_TWO_CMRES.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_TWO_CMRES;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_1024_A.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_1024_A;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_19937_A.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_19937_A;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_19937_C.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_19937_C;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_44497_A.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_44497_A;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_44497_B.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_44497_B;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_WELL_512_A.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_WELL_512_A;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.ACRNG_XOR_SHIFT_1024_S.toString())) {
                tmpRandomType = Factory.RandomType.ACRNG_XOR_SHIFT_1024_S;
            } else if (aRandomNumberGeneratorTypeRepresentation.equals(Factory.RandomType.PCG_32.toString())) {
                tmpRandomType = Factory.RandomType.PCG_32;
            }
            return tmpRandomType;
        }
        
    }  

    /**
     * DPD type
     */
    public enum DpdType {
        
        /**
         * DPD with cut-off length of 1.0
         */
        CUTOFF_LENGTH_ONE
        
    }
    
    /**
     * Electrostatics type
     */
    public enum ElectrostaticsType {
        
        /**
         * Ad-hoc electrostatics
         */
        AD_HOC
        
    }
    
    /**
     * Harmonic bond type
     */
    public enum HarmonicBondType {
        
        /**
         * Default harmonic bond
         */
        DEFAULT
        
    }

    /**
     * Integration type
     */
    public enum IntegrationType {
            
        /**
         * Groot-Warren Modified Velocity-Verlet (GWMVV)
         */
        GWMVV,
        /**
         * Self-consistent Modified Velocity-Verlet (SCMVV)
         */
        SCMVV,
        /**
         * Shardlow S1 Modified Velocity-Verlet (S1MVV)
         */
        S1MVV,
        /**
         * Nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN)
         */
        PNHLN;

        /**
         * Integration type  representations
         * 
         * @return Integration type representations
         */
        public static String[] getIntegrationTypeRepresentations() {
            return new String[] {
                IntegrationType.GWMVV.toString(),
                IntegrationType.SCMVV.toString(),
                IntegrationType.S1MVV.toString(),
                IntegrationType.PNHLN.toString()
            };
        }

        /**
         * Default integration type representation
         * 
         * @return Default integration type representation
         */
        public static String getDefaultIntegrationTypeRepresentation() {
            return Factory.IntegrationType.GWMVV.toString();
        }
        
    }    
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Random number generator type
     */
    private final RandomType randomType;

    /**
     * Number of random number generator warm-up steps
     */
    private final int numberORandomNumberGeneratorfWarmUpSteps;
    
    /**
     * DPD type
     */
    private final DpdType dpdType;

    /**
     * Electrostatics type
     */
    private final ElectrostaticsType electrostaticsType;
    
    /**
     * Harmonic bond type
     */
    private final HarmonicBondType harmonicBondType;
    
    /**
     * Integration type
     */
    private final IntegrationType integrationType;
    
    /**
     * Integration parameters
     */
    private final Object[] integrationParameters;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aRandomType Random number generator type
     * @param aNumberORandomNumberGeneratorfWarmUpSteps Number of random number 
     * generator warm-up steps (greater/equal 0)
     * @param aDpdType DPD type
     * @param anElectrostaticsType Electrostatics type
     * @param aHarmonicBondType Harmonic bond type
     * @param anIntegrationType Integration type
     * @param anIntegrationParameters Integration parameters
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public Factory(
        RandomType aRandomType,
        int aNumberORandomNumberGeneratorfWarmUpSteps,
        DpdType aDpdType,
        ElectrostaticsType anElectrostaticsType,
        HarmonicBondType aHarmonicBondType,
        IntegrationType anIntegrationType,
        Object[] anIntegrationParameters
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aRandomType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aRandomType is null.");
        }
        if (aNumberORandomNumberGeneratorfWarmUpSteps < 0) {
            throw new IllegalArgumentException("Factory.Constructor: aNumberORandomNumberGeneratorfWarmUpSteps < 0.");
        }
        if (aDpdType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aDpdType is null.");
        }
        if (anElectrostaticsType == null) {
            throw new IllegalArgumentException("Factory.Constructor: anElectrostaticsType is null.");
        }
        if (aHarmonicBondType == null) {
            throw new IllegalArgumentException("Factory.Constructor: aHarmonicBondType is null.");
        }
        if (anIntegrationType == null) {
            throw new IllegalArgumentException("Factory.Constructor: anIntegrationType is null.");
        }
        if (anIntegrationType == IntegrationType.GWMVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.GWMVV.");
            }
            if (anIntegrationParameters.length < 1) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 1 for IntegrationType.GWMVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Double)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Double for IntegrationType.GWMVV.");
            }
        }
        if (anIntegrationType == IntegrationType.SCMVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.SCMVV.");
            }
            if (anIntegrationParameters.length < 2) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 2 for IntegrationType.SCMVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Integer)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Integer for IntegrationType.SCMVV.");
            }
            if (!(anIntegrationParameters[1] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[1] is NOT instance of Boolean for IntegrationType.SCMVV.");
            }
        }
        if (anIntegrationType == IntegrationType.S1MVV) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.S1MVV.");
            }
            if (anIntegrationParameters.length < 1) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 1 for IntegrationType.S1MVV.");
            }
            if (!(anIntegrationParameters[0] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Boolean for IntegrationType.S1MVV.");
            }
        }
        if (anIntegrationType == IntegrationType.PNHLN) {
            if (anIntegrationParameters == null) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters is null for IntegrationType.PNHLN.");
            }
            if (anIntegrationParameters.length < 2) {
                throw new IllegalArgumentException("Factory.Constructor: Length of anIntegrationParameters < 2 for IntegrationType.PNHLN.");
            }
            if (!(anIntegrationParameters[0] instanceof Double)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[0] is NOT instance of Double for IntegrationType.PNHLN.");
            }
            if (!(anIntegrationParameters[1] instanceof Boolean)) {
                throw new IllegalArgumentException("Factory.Constructor: anIntegrationParameters[1] is NOT instance of Boolean for IntegrationType.PNHLN.");
            }
        }
        // </editor-fold>
        this.randomType = aRandomType;
        this.numberORandomNumberGeneratorfWarmUpSteps = aNumberORandomNumberGeneratorfWarmUpSteps;
        this.dpdType = aDpdType;
        this.electrostaticsType = anElectrostaticsType;
        this.harmonicBondType = aHarmonicBondType;
        this.integrationType = anIntegrationType;
        this.integrationParameters = anIntegrationParameters;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    // <editor-fold defaultstate="collapsed" desc="- Random number generators">
    /**
     * Returns new random number generator instance
     * 
     * @param aSeed Seed value (greater/equal 0)
     * @param this.numberORandomNumberGeneratorfWarmUpSteps Number of warm-up steps (greater/equal 0)
     * @return New random number generator instance
     */
    public IRandom getNewRandomNumberGenerator(int aSeed) {
        switch(this.randomType) {
            case ACRNG_ISAAC:
                return new ApacheCommonsRandom(RandomSource.ISAAC, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_JDK:
                return new ApacheCommonsRandom(RandomSource.JDK, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_KISS:
                return new ApacheCommonsRandom(RandomSource.KISS, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_MT:
                return new ApacheCommonsRandom(RandomSource.MT, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_MT_64:
                return new ApacheCommonsRandom(RandomSource.MT_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_MWC_256:
                return new ApacheCommonsRandom(RandomSource.MWC_256, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_SPLIT_MIX_64:
                return new ApacheCommonsRandom(RandomSource.SPLIT_MIX_64, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_TWO_CMRES:
                return new ApacheCommonsRandom(RandomSource.TWO_CMRES, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_1024_A:
                return new ApacheCommonsRandom(RandomSource.WELL_1024_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_19937_A:
                return new ApacheCommonsRandom(RandomSource.WELL_19937_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_19937_C:
                return new ApacheCommonsRandom(RandomSource.WELL_19937_C, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_44497_A:
                return new ApacheCommonsRandom(RandomSource.WELL_44497_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_44497_B:
                return new ApacheCommonsRandom(RandomSource.WELL_44497_B, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_WELL_512_A:
                return new ApacheCommonsRandom(RandomSource.WELL_512_A, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case ACRNG_XOR_SHIFT_1024_S:
                return new ApacheCommonsRandom(RandomSource.XOR_SHIFT_1024_S, aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            case PCG_32:
                return new PcgRandom(aSeed, this.numberORandomNumberGeneratorfWarmUpSteps);
            default:
                throw new IllegalArgumentException("Factory.getNewRandomNumberGenerator: Unknown random type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- DPD calculators">
    // <editor-fold defaultstate="collapsed" desc="-- PNHLN calculators">
    /**
     * Returns particle pair velocity update plus G calculator for PNHLN integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair velocity update plus G calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionPnhlnCalculator getParticlePairPnhlnVelocityUpdatePlusGCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairPnhlnVelocityUpdatePlusGCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair velocity update plus G calculator for PNHLN integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair velocity update plus G calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionPnhlnCalculator getParticlePairPnhlnVelocityUpdatePlusGCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairPnhlnVelocityUpdatePlusGCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairPnhlnVelocityUpdatePlusGCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- S1MVV calculators">
    /**
     * Returns particle pair velocity update calculator for S1MVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair velocity update calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairS1mvvVelocityUpdateCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairS1mvvVelocityUpdateCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairS1mvvVelocityUpdateCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair velocity update calculator for S1MVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair velocity update calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairS1mvvVelocityUpdateCalculator(
        IParticlePairInteractionCalculator aParticlePairInteractionCalculator
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairS1mvvVelocityUpdateCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairS1mvvVelocityUpdateCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- SCMVV calculators">
    /**
     * Returns particle pair DPD full force calculator for SCMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceFullCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceFullCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD full force calculator for SCMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceFullCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceFullCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD dissipative force calculator for 
     * SCMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceDissipativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD dissipative force calculator for 
     * SCMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairScmvvDpdForceDissipativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairScmvvDpdForceDissipativeCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairScmvvDpdForceDissipativeCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- GWMVV calculators">
    /**
     * Returns particle pair DPD full force calculator for GWMVV integration type
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairGwmvvDpdForceFullCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairGwmvvDpdForceFullCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairGwmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD full force calculator for GWMVV integration type
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairGwmvvDpdForceFullCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairGwmvvDpdForceFullCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairGwmvvDpdForceFullCalculator: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- General calculators">
    /**
     * Returns particle pair DPD conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceConservativeCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceConservativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD conservative force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD conservative force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairDpdForceConservativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdForceConservativeCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdForceConservativeCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairDpdPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdPotentialCutoff1Calculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdPotentialCalculator: Unknown DPD type.");
        }
    }

    /**
     * Returns particle pair DPD potential calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairDpdPotentialCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                return new ParticlePairDpdPotentialCutoff1Calculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairDpdPotentialCalculator: Unknown DPD type.");
        }
    }


    /**
     * Returns particle-pair nearest-neighbor calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair DPD potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairNearestNeighborCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException 
    {
        return new ParticlePairNearestNeighborCalculator(
            this,
            aSimulationLogger, 
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength,
            aParallelizationInfo,
            aRandomNumberSeed
        );
    }
    // </editor-fold>
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- DPD parameters">
    /**
     * DPD cut-off length (in DPD units)
     * 
     * @return DPD cut-off length (in DPD units)
     */
    public double getDpdCutOffLength() {
        switch (this.dpdType) {
            case CUTOFF_LENGTH_ONE:
                // DPD cut-off length of 1.0
                return 1.0;
            default:
                throw new IllegalArgumentException("Factory.getDpdCutOffLength: Unknown DPD type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Harmonic bond calculators">
    /**
     * Returns harmonic bond conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @return Harmonic bond force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondForceCalculator getHarmonicBondForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        switch (this.harmonicBondType) {
            case DEFAULT:
                return new HarmonicBondForceConservativeCalculator(
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aParallelizationInfo
                );
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondForceConservativeCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond conservative force calculator
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @return Harmonic bond force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondForceCalculator getHarmonicBondForceConservativeCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        switch (this.harmonicBondType) {
            case DEFAULT:
                return new HarmonicBondForceConservativeCalculator(aHarmonicBondPropertyCalculator);
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondForceConservativeCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @return Harmonic bond potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondPropertyCalculator getHarmonicBondPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        switch (this.harmonicBondType) {
            case DEFAULT:
                return new HarmonicBondPotentialCalculator(
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aParallelizationInfo
                );
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondPotentialCalculator: Unknown harmonic bond type.");
        }
    }

    /**
     * Returns harmonic bond potential calculator
     * 
     * @param aHarmonicBondPropertyCalculator HarmonicBondPropertyCalculator instance
     * @return Harmonic bond potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IHarmonicBondPropertyCalculator getHarmonicBondPotentialCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        switch (this.harmonicBondType) {
            case DEFAULT:
                return new HarmonicBondPotentialCalculator(aHarmonicBondPropertyCalculator);
            default:
                throw new IllegalArgumentException("Factory.getHarmonicBondPotentialCalculator: Unknown harmonic bond type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Electrostatics calculators">
    /**
     * Returns particle pair electrostatics conservative force calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair electrostatics force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairElectrostaticsForceConservativeCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocForceConservativeCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsForceConservativeCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics conservative force calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair electrostatics force calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairForceCalculator getParticlePairElectrostaticsForceConservativeCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocForceConservativeCalculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsForceConservativeCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics potential calculator
     * 
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @return Particle pair electrostatics potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairElectrostaticsPotentialCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocPotentialCalculator(
                    this,
                    aSimulationLogger, 
                    aBoxSize, 
                    aPeriodicBoundaries, 
                    aCutOffLength,
                    aParallelizationInfo,
                    aRandomNumberSeed
                );
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsPotentialCalculator: Unknown electrostatics type.");
        }
    }

    /**
     * Returns particle pair electrostatics potential calculator
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @return Particle pair electrostatics potential calculator
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public IParticlePairInteractionCalculator getParticlePairElectrostaticsPotentialCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        switch (this.electrostaticsType) {
            case AD_HOC:
                return new ParticlePairElectrostaticsAdHocPotentialCalculator(aParticlePairInteractionCalculator);
            default:
                throw new IllegalArgumentException("Factory.getParticlePairElectrostaticsPotentialCalculator: Unknown electrostatics type.");
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Time step calculators">
    /**
     * Time step calculator
     * 
     * @param aSimulationOutput Simulation output
     * @param aSimulationLogger Simulation logger
     * @param aParameters Parameters instance
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @param aParticlePositionPool Particle position pool (NOT allowed to be 
     * null)
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     * @return Time step calculator
     * @throws IllegalArgumentException Thrown if argument is invalid
     */
    public ITimeStepCalculator getTimeStepCalculator(
        IOutput aSimulationOutput,
        ILogger aSimulationLogger, 
        Parameters aParameters, 
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed,
        ParticlePositionPool aParticlePositionPool,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        boolean tmpIsCacheUsage = false;
        switch (this.integrationType) {
            case GWMVV:
                // Lambda parameter for Groot-Warren Modified Velocity-Verlet (GWMVV) integration
                double tmpLambda = (Double) this.integrationParameters[0];
                return new GwmvvTimeStepCalculator(
                    this,
                    tmpLambda,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case SCMVV:
                // Number of self-constistent iterations for Self-consistent Modified Velocity-Verlet (SCMVV) integration
                int tmpSelfConsistentIterationNumber = (Integer) this.integrationParameters[0];
                // Flag for use of cache for Self-consistent Modified Velocity-Verlet (SCMVV) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[1];
                return new ScmvvTimeStepCalculator(
                    this,
                    tmpSelfConsistentIterationNumber,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case S1MVV:
                // Flag for use of cache for Shardlow S1 Modified Velocity-Verlet (S1MVV) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[0];
                return new S1mvvTimeStepCalculator(
                    this,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            case PNHLN:
                // Mu parameter for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
                double tmpMu = (Double) this.integrationParameters[0];
                // Flag for use of cache for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
                tmpIsCacheUsage = (Boolean) this.integrationParameters[1];
                return new PnhlnTimeStepCalculator(
                    this,
                    tmpMu,
                    tmpIsCacheUsage,
                    aSimulationOutput,
                    aSimulationLogger, 
                    aParameters, 
                    aParallelizationInfo,
                    aRandomNumberSeed,
                    aParticlePositionPool,
                    aMaximumNumberOfPositionCorrectionTrials
                );
            default:
                throw new IllegalArgumentException("Factory.getTimeStepCalculator: Unknown integration type.");
        }
    }
    // </editor-fold>
    // </editor-fold>
    
}
