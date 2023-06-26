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
package de.gnwi.jdpdsp.samples;

import de.gnwi.jdpdsp.utilities.Factory;
import de.gnwi.jdpdsp.samples.harmonicBonds.ParticlePairHarmonicBond;
import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBond;
import java.io.File;
import java.util.HashMap;
import de.gnwi.jdpdsp.interfaces.IInput;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.movement.MoleculeAccelerationDescription;
import de.gnwi.jdpdsp.movement.MoleculeBoundaryDescription;
import de.gnwi.jdpdsp.parameters.ParticleTypes;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.LinkedList;
import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.utilities.Constants;
import de.gnwi.jdpdsp.utilities.Electrostatics;
import de.gnwi.jdpdsp.utilities.MoleculeDescription;
import de.gnwi.jdpdsp.movement.MoleculeFixationDescription;
import de.gnwi.jdpdsp.movement.MoleculeSphereDescription;
import de.gnwi.jdpdsp.rg.MoleculeRgCalculationDescription;
import de.gnwi.jdpdsp.movement.MoleculeVelocityFixationDescription;
import de.gnwi.jdpdsp.nearestNeighbor.NearestNeighborBaseParticleDescription;
import de.gnwi.jdpdsp.utilities.GravitationalAcceleration;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import de.gnwi.jdpdsp.utilities.Utils;
import java.util.TreeMap;

/**
 * File input
 * @author Achim Zielesny
 */
public class FileInput implements IInput {

    // <editor-fold defaultstate="collapsed" desc="Enum Section">
    /**
     * Log level
     */
    private enum Section {
            
        /**
         * Section GENERAL
         */
        GENERAL,
        /**
         * Section PARTICLE_DESCRIPTION
         */
        PARTICLE_DESCRIPTION,
        /**
         * Section CHEMICAL_SYSTEM_DESCRIPTION
         */
        CHEMICAL_SYSTEM_DESCRIPTION,
        /**
         * Section INTERACTION_DESCRIPTION
         */
        INTERACTION_DESCRIPTION,
        /**
         * Section SIMULATION_DESCRIPTION
         */
        SIMULATION_DESCRIPTION,
        /**
         * Section SIMULATION_COUNTS
         */
        SIMULATION_COUNTS

    }    
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    // <editor-fold defaultstate="collapsed" desc="- Comment related strings">
    /**
     * Comment line prefix
     */
    private static final String COMMENT_LINE_PREFIX = "#";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Section related strings">
    /**
     * Format string for section start tag
     */
    private static final String SECTION_START_TAG_FORMAT = "[%s]";

    /**
     * Format string for section end tag
     */
    private static final String SECTION_END_TAG_FORMAT = "[/%s]";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Table related strings">
    /**
     * Table start string
     */
    private static final String TABLE_START = "TABLE_START";

    /**
     * Table end string
     */
    private static final String TABLE_END = "TABLE_END";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Properties related strings">
    // <editor-fold defaultstate="collapsed" desc="-- General properties">
    /**
     * String for parameter Version
     */
    private static final String VERSION = "Version";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- Particle description properties">
    /**
     * String for parameter ParticleTable
     */
    private final static String PARTICLE_TABLE = "ParticleTable";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- Chemical system description properties">
    /**
     * String for parameter MoleculeTable
     */
    private static final String MOLECULE_TABLE = "MoleculeTable";

    /**
     * String for parameter BoxSize
     */
    private static final String BOX_SIZE = "BoxSize";
    
    /**
     * String for parameter MoleculeFixation
     */
    private static final String MOLECULE_FIXATION = "MoleculeFixation";
    
    /**
     * String for parameter MoleculeBoundary
     */
    private static final String MOLECULE_BOUNDARY = "MoleculeBoundary";
    
    /**
     * String for parameter MoleculeSphere
     */
    private static final String MOLECULE_SPHERE = "MoleculeSphere";

    /**
     * String for parameter MoleculeFixedVelocity
     */
    private static final String MOLECULE_FIXED_VELOCITY = "MoleculeFixedVelocity";

    /**
     * String for parameter MoleculeAcceleration
     */
    private static final String MOLECULE_ACCELERATION = "MoleculeAcceleration";

    /**
     * String for parameter RadiusOfGyration
     */
    private static final String RADIUS_OF_GYRATION = "RadiusOfGyration";

    /**
     * String for parameter NearestNeighborParticle
     */
    private static final String NEAREST_NEIGHBOR_PARTICLE = "NearestNeighborParticle";
    
    /**
     * String for parameter NearestNeighborDistance
     */
    private static final String NEAREST_NEIGHBOR_DISTANCE = "NearestNeighborDistance";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- Interaction description properties">
    /**
     * String for parameter Temperature
     */
    private static final String TEMPERATURE = "Temperature";
    
    /**
     * String for parameter DpdSigma
     */
    private static final String DPD_SIGMA = "DpdSigma";
    
    /**
     * String for parameter InteractionTable
     */
    private static final String INTERACTION_TABLE = "InteractionTable";
    
    /**
     * String for parameter Bonds12Table
     */
    private static final String BONDS_12_TABLE = "Bonds12Table";
    
    /**
     * String for parameter IsGaussianRandomDpdForce
     */
    private static final String IS_GAUSSIAN_RANDOM_DPD_FORCE = "IsGaussianRandomDpdForce";
    
    /**
     * String for parameter Electrostatics
     */
    private static final String ELECTROSTATICS = "Electrostatics";
    
    /**
     * String for parameter GravitationalAcceleration
     */
    private static final String GRAVITATIONAL_ACCELERATION = "GravitationalAcceleration";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- Simulation description properties">
    /**
     * String for parameter TimeStepNumber
     */
    private static final String TIME_STEP_NUMBER = "TimeStepNumber";

    /**
     * String for parameter TimeStepLength
     */
    private static final String TIME_STEP_LENGTH = "TimeStepLength";

    /**
     * String for parameter TimeStepFrequencyForOutput
     */
    private static final String TIME_STEP_FREQUENCY_FOR_OUTPUT = "TimeStepFrequencyForOutput";

    /**
     * String for parameter IntegrationType
     */
    private static final String INTEGRATION_TYPE = "IntegrationType";

    /**
     * String for parameter InitialPotentialEnergyMinimizationStepNumber
     */
    private static final String INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_NUMBER = "InitialPotentialEnergyMinimizationStepNumber";
    
    /**
     * String for parameter IsInitialPotentialEnergyMinimizationStepOutput
     */
    private static final String IS_INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_OUTPUT = "IsInitialPotentialEnergyMinimizationStepOutput";
    
    /**
     * String for parameter PeriodicBoundaries
     */
    private static final String PERIODIC_BOUNDARIES = "PeriodicBoundaries";

    /**
     * String for parameter IsDpdUnitMass
     */
    private static final String IS_DPD_UNIT_MASS = "IsDpdUnitMass";

    /**
     * String for parameter IsVelocityScaling
     */
    private static final String IS_VELOCITY_SCALING = "IsVelocityScaling";

    /**
     * String for parameter InitialVelocityScalingSteps
     */
    private static final String INITIAL_VELOCITY_SCALING_STEPS = "InitialVelocityScalingSteps";

    /**
     * String for parameter RandomNumberGenerator
     */
    private static final String RANDOM_NUMBER_GENERATOR = "RandomNumberGenerator";
    
    /**
     * Default random number generator warm-up steps
     */
    private static final int DEFAULT_RANDOM_WARM_UP_STEPS = 10000;
    
    /**
     * Default random seed
     */
    private static final int DEFAULT_RANDOM_SEED = 123456;
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="-- Simulation counts properties">
    /**
     * Particle number
     */
    private static final String PARTICLE_NUMBER = "ParticleNumber";
    // </editor-fold>
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    // <editor-fold defaultstate="collapsed" desc="- Input file related variables">
    /**
     * Pathname of input file
     */
    private final String inputFilePathname;
    
    /**
     * Full path of input
     */
    private final String inputPath;
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Simulation logger">
    /**
     * Simulation logger
     */
    private ILogger simulationLogger;
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Section to content map">
    /**
     * Map that maps section tag to its jagged string array content
     */
    private final HashMap<Section, String[][]> sectionMap;
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Particle token map">
    private final TreeMap<String, String> particleTokenMap;
    // </editor-fold>
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
     * @param anInputFilePathname Pathname of input file
     * @param aSimulationLogger Simulation logger
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public FileInput(String anInputFilePathname, ILogger aSimulationLogger) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anInputFilePathname == null || anInputFilePathname.isEmpty() || !(new File(anInputFilePathname)).isFile()) {
            throw new IllegalArgumentException("FileInput.Constructor: anInputFilePathname is null/empty or does not exist.");
        }
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("FileInput.Constructor: aSimulationLogger is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        // Default value for this.maximumNumberOfPositionCorrectionTrials is arbitrarily 
        // chosen to be 1000 for all practical purposes
        this.maximumNumberOfPositionCorrectionTrials = 1000;
        this.particleTokenMap = new TreeMap<>();
        this.inputFilePathname = anInputFilePathname;
        try {
            this.inputPath = (new File(anInputFilePathname)).getParent();
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.Constructor: Can not set input path.");
        }
        this.sectionMap = new HashMap<>(Section.values().length);
        if (!this.readInputFile()) {
            throw new IllegalArgumentException("FileInput.Constructor: Can not read input file.");
        }
        // Check version AFTER input file is read
        String tmpVersion = this.getSingleStringValue(Section.GENERAL, FileInput.VERSION);
        if (tmpVersion == null || tmpVersion.isEmpty() || !tmpVersion.equals("1.0.0.0")) {
            throw new IllegalArgumentException("FileInput.Constructor: Illegal version of input file.");
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    // <editor-fold defaultstate="collapsed" desc="- Factory (get)">
    /**
     * Returns factory
     * 
     * @return Factory
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    @Override
    public Factory getFactory() {
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = this.simulationLogger.getId();
        this.simulationLogger.appendMethodCallStart("FileInput.getFactory", tmpId);
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Random number generator">
        Factory.RandomType tmpRandomType = Factory.RandomType.getDefaultRandomNumberGeneratorType();
        int tmpRandomSeed = FileInput.DEFAULT_RANDOM_SEED;
        int tmpNumberORandomNumberGeneratorWarmUpSteps = FileInput.DEFAULT_RANDOM_WARM_UP_STEPS;
        if (this.hasParameter(Section.SIMULATION_DESCRIPTION, FileInput.RANDOM_NUMBER_GENERATOR)) {
            String[] tmpValues = this.getTripleStringValues(Section.SIMULATION_DESCRIPTION, FileInput.RANDOM_NUMBER_GENERATOR);
            // Check for legacy definitions: Replace invalid old settings by default
            if (Factory.RandomType.isDefinedRandomNumberGeneratorTypeRepresentation(tmpValues[0])) {
                tmpRandomType = Factory.RandomType.valueOf(tmpValues[0]);
            }
            tmpRandomSeed = Integer.valueOf(tmpValues[1]);
            tmpNumberORandomNumberGeneratorWarmUpSteps = Integer.valueOf(tmpValues[2]);
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Integration type">
        Factory.IntegrationType tmpIntegrationType = null;
        Object[] tmpIntegrationParameters = null;
        String[] tmpValues = this.getMultipleStringValues(Section.SIMULATION_DESCRIPTION, FileInput.INTEGRATION_TYPE);
        if (tmpValues[0].equals(Factory.IntegrationType.GWMVV.toString())) {
            tmpIntegrationType = Factory.IntegrationType.GWMVV;
            // Float parameter: Lambda parameter for Groot-Warren Modified Velocity-Verlet (GWMVV) integration
            tmpIntegrationParameters = 
                new Object[] {
                    Float.valueOf(tmpValues[1])
                };
        } else if (tmpValues[0].equals(Factory.IntegrationType.SCMVV.toString())) {
            tmpIntegrationType = Factory.IntegrationType.SCMVV;
            // Integer parameter: Number of self-constistent iterations for Self-consistent Modified Velocity-Verlet (SCMVV) integration
            // Boolean parameter: Flag for use of cache for Self-consistent Modified Velocity-Verlet (SCMVV) integration
            tmpIntegrationParameters = 
                new Object[] {
                    Integer.valueOf(tmpValues[1]),
                    Boolean.valueOf(tmpValues[2])
                };
        } else if (tmpValues[0].equals(Factory.IntegrationType.S1MVV.toString())) {
            tmpIntegrationType = Factory.IntegrationType.S1MVV;
            // Boolean parameter: Flag for use of cache for Shardlow S1 Modified Velocity-Verlet (S1MVV) integration
            tmpIntegrationParameters = 
                new Object[] {
                    Boolean.valueOf(tmpValues[1])
                };
        } else if (tmpValues[0].equals(Factory.IntegrationType.PNHLN.toString())) {
            tmpIntegrationType = Factory.IntegrationType.PNHLN;
            // Float parameter: Mu parameter for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
            // Boolean parameter: Flag for use of cache for nonsymmetric pairwise Nose–Hoover–Langevin thermostat (PNHLN) integration
            tmpIntegrationParameters = 
                new Object[] {
                    Float.valueOf(tmpValues[1]),
                    Boolean.valueOf(tmpValues[2])
                };
        } else {
            throw new IllegalArgumentException("FileInput.getFactory: Illegal integration type");
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        Factory.ElectrostaticsType tmpElectrostaticsType = Factory.ElectrostaticsType.AD_HOC;
        if (this.isElectrostaticsTypeDpd()) {
            tmpElectrostaticsType = Factory.ElectrostaticsType.DPD;
        }
        // </editor-fold>
        Factory tmpFactory = 
            new Factory(
                tmpRandomType, 
                tmpNumberORandomNumberGeneratorWarmUpSteps,
                Factory.DpdType.CUTOFF_LENGTH_ONE, 
                tmpElectrostaticsType,
                Factory.BondType.HARMONIC,
                tmpIntegrationType,
                tmpIntegrationParameters
            );
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Random type           = " + tmpRandomType.toString());
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Random seed           = " + String.valueOf(tmpRandomSeed));
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Random warm-up steps  = " + String.valueOf(tmpNumberORandomNumberGeneratorWarmUpSteps));
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: DPD type              = " + Factory.DpdType.CUTOFF_LENGTH_ONE.toString());
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Electrostatics type   = " + tmpElectrostaticsType.toString());
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Harmonic bond type    = " + Factory.BondType.HARMONIC.toString());
        this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Integration type      = " + tmpIntegrationType.toString());
        if (tmpIntegrationParameters != null) {
            for (Object tmpParameter : tmpIntegrationParameters) {
                this.simulationLogger.appendIntermediateResults("FileInput.getFactory: Integration parameter = " + tmpParameter.toString());
            }
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCallEnd("FileInput.getFactory", tmpId);
        // </editor-fold>
        return tmpFactory;
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Particle description properties (get)">
    /**
     * Particle types
     * 
     * @return Particle types
     */
    @Override
    public ParticleTypes getParticleTypes() {
        String[][] tmpParticleTable = this.getValueTable(Section.PARTICLE_DESCRIPTION, FileInput.PARTICLE_TABLE);
        HashMap<String, Integer> tmpParticleTokenToIndexMap = new HashMap<>(tmpParticleTable.length);
        String[] aParticleTokens = new String[tmpParticleTable.length];
        float[] aMolarMasses = new float[tmpParticleTable.length];
        float[] aCharges = new float[tmpParticleTable.length];
        for (int i = 0; i < tmpParticleTable.length; i++) {
            tmpParticleTokenToIndexMap.put(tmpParticleTable[i][0], i);
            aParticleTokens[i] = tmpParticleTable[i][0];
            aCharges[i] = Float.valueOf(tmpParticleTable[i][1]);
            aMolarMasses[i] = Float.valueOf(tmpParticleTable[i][2]);
        }
        return new ParticleTypes(
            tmpParticleTokenToIndexMap, 
            aParticleTokens, 
            aCharges, 
            aMolarMasses
        );
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Chemical system description properties (get)">
    /**
     * Molecule descriptions
     * 
     * @return Molecule descriptions
     */
    @Override
    public MoleculeDescription[] getMoleculeDescriptions() {
        String[][] tmpMoleculeTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_TABLE);
        MoleculeDescription[] tmpMoleculeDescriptions = new MoleculeDescription[tmpMoleculeTable.length];
        for (int i = 0; i < tmpMoleculeTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeTable[i][0];
            String tmpPositionsBondsFilePathname = this.inputPath + File.separatorChar + tmpMoleculeTable[i][1];
            String[][] tmpPositionsBondsFileContent = this.readPositionsBondsFile(tmpPositionsBondsFilePathname);
            int tmpTotalMoleculeParticleNumber = Integer.valueOf(tmpPositionsBondsFileContent[0][1]);
            int tmpSingleMoleculeParticleNumber = Integer.valueOf(tmpPositionsBondsFileContent[1][1]);
            String[] tmpParticleTokens = new String[tmpTotalMoleculeParticleNumber];
            int[] tmpParticleBackboneIndices = new int[tmpTotalMoleculeParticleNumber];
            float[] tmpR_x = new float[tmpTotalMoleculeParticleNumber];
            float[] tmpR_y = new float[tmpTotalMoleculeParticleNumber];
            float[] tmpR_z = new float[tmpTotalMoleculeParticleNumber];
            int[][] tmpBondOffsets = new int[tmpTotalMoleculeParticleNumber][];
            for (int j = 0; j < tmpTotalMoleculeParticleNumber; j++) {
                // Use referenced particle token to avoid multiple storage of same particle token
                tmpParticleTokens[j] = this.getReferencedParticleToken(tmpPositionsBondsFileContent[j + 2][1]);
                tmpParticleBackboneIndices[j] = Integer.valueOf(tmpPositionsBondsFileContent[j + 2][2]);
                tmpR_x[j] = Float.valueOf(tmpPositionsBondsFileContent[j + 2][3]);
                tmpR_y[j] = Float.valueOf(tmpPositionsBondsFileContent[j + 2][4]);
                tmpR_z[j] = Float.valueOf(tmpPositionsBondsFileContent[j + 2][5]);
                // Read bonds
                if (tmpPositionsBondsFileContent[j + 2].length > 6) {
                    int[] tmpSingleParticleBondOffsets = new int[tmpPositionsBondsFileContent[j + 2].length - 6];
                    for (int k = 0; k < tmpSingleParticleBondOffsets.length; k++) {
                        tmpSingleParticleBondOffsets[k] = Integer.valueOf(tmpPositionsBondsFileContent[j + 2][k + 6]);
                    }
                    tmpBondOffsets[j] = tmpSingleParticleBondOffsets;
                } else {
                    tmpBondOffsets[j] = null;
                }
            }
            // Read backbone bonds if available
            HarmonicBond[] tmpBackboneBonds = null;
            if (tmpPositionsBondsFileContent.length > tmpTotalMoleculeParticleNumber + 2) {
                int tmpBackboneBondNumber = Integer.valueOf(tmpPositionsBondsFileContent[tmpTotalMoleculeParticleNumber + 2][1]);
                tmpBackboneBonds = new HarmonicBond[tmpBackboneBondNumber];
                int tmpOffset = tmpTotalMoleculeParticleNumber + 3;
                for (int j = 0; j < tmpBackboneBondNumber; j++) {
                    if (tmpPositionsBondsFileContent[tmpOffset + j].length == 4) {
                        // NOTE: Use HarmonicBond.HarmonicBondBehaviour.DEFAULT for backbone bonds
                        //       without HarmonicBondBehaviour specification
                        tmpBackboneBonds[j] = 
                            new HarmonicBond(
                                Integer.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][0]),
                                Integer.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][1]),
                                Float.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][2]),
                                Float.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][3]),
                                HarmonicBond.HarmonicBondBehaviour.DEFAULT
                            );
                    } else {
                        tmpBackboneBonds[j] = 
                            new HarmonicBond(
                                Integer.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][0]),
                                Integer.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][1]),
                                Float.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][2]),
                                Float.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][3]),
                                HarmonicBond.HarmonicBondBehaviour.valueOf(tmpPositionsBondsFileContent[tmpOffset + j][4])
                            );
                    }
                }
            }
            // Set molecule description
            tmpMoleculeDescriptions[i] = 
                new MoleculeDescription(
                    tmpMoleculeName,
                    tmpTotalMoleculeParticleNumber,
                    tmpSingleMoleculeParticleNumber,
                    tmpParticleTokens,
                    tmpParticleBackboneIndices,
                    tmpR_x,
                    tmpR_y,
                    tmpR_z,
                    tmpBondOffsets,
                    tmpBackboneBonds
                );
        }
        return tmpMoleculeDescriptions;
    }
    
    /**
     * Box size in DPD units
     * 
     * @return Box size in DPD units
     */
    @Override
    public BoxSize getBoxSize() {
        float[] tmpBoxSizeParameters =  this.getTripleFloatValues(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.BOX_SIZE);
        return new BoxSize(tmpBoxSizeParameters[0], tmpBoxSizeParameters[1], tmpBoxSizeParameters[2]);
    }
    
    /**
     * Molecule fixation descriptions
     * 
     * @return Molecule fixation descriptions or null if no molecule is fixed
     */
    @Override
    public MoleculeFixationDescription[] getMoleculeFixationDescriptions() {
        String[][] tmpMoleculeFixationTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_FIXATION);
        LinkedList<MoleculeFixationDescription> tmpMoleculeFixationList = new LinkedList<>();
        for (int i = 0; i < tmpMoleculeFixationTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeFixationTable[i][0];
            boolean tmpIsFixedX = Boolean.valueOf(tmpMoleculeFixationTable[i][1]);
            boolean tmpIsFixedY = Boolean.valueOf(tmpMoleculeFixationTable[i][2]);
            boolean tmpIsFixedZ = Boolean.valueOf(tmpMoleculeFixationTable[i][3]);

            // If tmpMaxTimeStep is NOT defined it gets the highest possible value
            int tmpMaxTimeStep = Constants.MAXIMUM_NUMBER_OF_TIME_STEPS;
            if (tmpMoleculeFixationTable[i].length > 4) {
                tmpMaxTimeStep = Integer.valueOf(tmpMoleculeFixationTable[i][4]);
            }
            
            if ((tmpIsFixedX || tmpIsFixedY || tmpIsFixedZ) && tmpMaxTimeStep > 0) {
                tmpMoleculeFixationList.add(
                    new MoleculeFixationDescription(
                        tmpMoleculeName, 
                        tmpIsFixedX, 
                        tmpIsFixedY, 
                        tmpIsFixedZ, 
                        tmpMaxTimeStep
                    )
                );
            }
        }
        if (!tmpMoleculeFixationList.isEmpty()) {
            return tmpMoleculeFixationList.toArray(new MoleculeFixationDescription[0]);
        } else {
            return null;
        }
    }
    
    /**
     * Molecule boundary descriptions
     * 
     * @return Molecule boundary descriptions or null if no molecule is bounded
     */
    @Override
    public MoleculeBoundaryDescription[] getMoleculeBoundaryDescriptions() {
        String[][] tmpMoleculeBoundaryTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_BOUNDARY);
        LinkedList<MoleculeBoundaryDescription> tmpMoleculeBoundaryList = new LinkedList<>();
        for (int i = 0; i < tmpMoleculeBoundaryTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeBoundaryTable[i][0];
            boolean tmpIsActiveX = Boolean.valueOf(tmpMoleculeBoundaryTable[i][1]);
            boolean tmpIsActiveY = Boolean.valueOf(tmpMoleculeBoundaryTable[i][4]);
            boolean tmpIsActiveZ = Boolean.valueOf(tmpMoleculeBoundaryTable[i][7]);

            // If tmpMaxTimeStep is NOT defined it gets the highest possible value
            int tmpMaxTimeStep = Constants.MAXIMUM_NUMBER_OF_TIME_STEPS;
            if (tmpMoleculeBoundaryTable[i].length > 10) {
                tmpMaxTimeStep = Integer.valueOf(tmpMoleculeBoundaryTable[i][10]);
            }
            
            if ((tmpIsActiveX || tmpIsActiveY || tmpIsActiveZ) && tmpMaxTimeStep > 0) {
                float tmpXmin = Float.valueOf(tmpMoleculeBoundaryTable[i][2]);
                float tmpXmax = Float.valueOf(tmpMoleculeBoundaryTable[i][3]);
                float tmpYmin = Float.valueOf(tmpMoleculeBoundaryTable[i][5]);
                float tmpYmax = Float.valueOf(tmpMoleculeBoundaryTable[i][6]);
                float tmpZmin = Float.valueOf(tmpMoleculeBoundaryTable[i][8]);
                float tmpZmax = Float.valueOf(tmpMoleculeBoundaryTable[i][9]);
                tmpMoleculeBoundaryList.add(
                    new MoleculeBoundaryDescription(
                        tmpMoleculeName, 
                        tmpIsActiveX,
                        tmpXmin,
                        tmpXmax,
                        tmpIsActiveY,
                        tmpYmin,
                        tmpYmax,
                        tmpIsActiveZ,
                        tmpZmin,
                        tmpZmax,
                        tmpMaxTimeStep
                    )
                );
            }
        }
        if (!tmpMoleculeBoundaryList.isEmpty()) {
            return tmpMoleculeBoundaryList.toArray(new MoleculeBoundaryDescription[0]);
        } else {
            return null;
        }
    }

    /**
     * Molecule sphere descriptions
     * 
     * @return Molecule sphere descriptions or null if no molecule has a exclusion/inclusion sphere
     */
    @Override
    public MoleculeSphereDescription[] getMoleculeSphereDescriptions() {
        if (this.hasParameter(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_SPHERE)) {
            String[][] tmpMoleculeSphereTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_SPHERE);
            LinkedList<MoleculeSphereDescription> tmpMoleculeSphereList = new LinkedList<>();
            for (int i = 0; i < tmpMoleculeSphereTable.length; i++) {
                float tmpSphereRadius = Float.valueOf(tmpMoleculeSphereTable[i][5]);
                if (tmpSphereRadius > 0.0f) {
                    String tmpMoleculeName = tmpMoleculeSphereTable[i][0];
                    boolean tmpIsExclusiveSphere = Boolean.valueOf(tmpMoleculeSphereTable[i][1]);
                    float tmpSphereCenterX = Float.valueOf(tmpMoleculeSphereTable[i][2]);
                    float tmpSphereCenterY = Float.valueOf(tmpMoleculeSphereTable[i][3]);
                    float tmpSphereCenterZ = Float.valueOf(tmpMoleculeSphereTable[i][4]);
                    int tmpMaxTimeStep = Integer.valueOf(tmpMoleculeSphereTable[i][6]);
                    tmpMoleculeSphereList.add(
                        new MoleculeSphereDescription(
                            tmpMoleculeName, 
                            tmpIsExclusiveSphere,
                            tmpSphereCenterX,
                            tmpSphereCenterY,
                            tmpSphereCenterZ,
                            tmpSphereRadius,
                            tmpMaxTimeStep
                        )
                    );
                }
            }
            if (!tmpMoleculeSphereList.isEmpty()) {
                return tmpMoleculeSphereList.toArray(new MoleculeSphereDescription[0]);
            } else {
                return null;
            }
        } else {
            return null;
        }
    }
    
    /**
     * Molecule velocity fixation descriptions
     * 
     * @return Molecule velocity fixation descriptions or null if no molecule velocity is fixed
     */
    @Override
    public MoleculeVelocityFixationDescription[] getMoleculeVelocityFixationDescriptions() {
        String[][] tmpMoleculeVelocityFixationTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_FIXED_VELOCITY);
        LinkedList<MoleculeVelocityFixationDescription> tmpMoleculeVelocityFixationList = new LinkedList<>();
        for (int i = 0; i < tmpMoleculeVelocityFixationTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeVelocityFixationTable[i][0];
            boolean tmpIsFixedX = Boolean.valueOf(tmpMoleculeVelocityFixationTable[i][1]);
            float tmpVelocityX = Float.valueOf(tmpMoleculeVelocityFixationTable[i][2]);
            boolean tmpIsFixedY = Boolean.valueOf(tmpMoleculeVelocityFixationTable[i][3]);
            float tmpVelocityY = Float.valueOf(tmpMoleculeVelocityFixationTable[i][4]);
            boolean tmpIsFixedZ = Boolean.valueOf(tmpMoleculeVelocityFixationTable[i][5]);
            float tmpVelocityZ = Float.valueOf(tmpMoleculeVelocityFixationTable[i][6]);

            // If tmpMaxTimeStep is NOT defined it gets the highest possible value
            int tmpMaxTimeStep = Constants.MAXIMUM_NUMBER_OF_TIME_STEPS;
            if (tmpMoleculeVelocityFixationTable[i].length > 7) {
                tmpMaxTimeStep = Integer.valueOf(tmpMoleculeVelocityFixationTable[i][7]);
            }
            
            if ((tmpIsFixedX || tmpIsFixedY || tmpIsFixedZ) && tmpMaxTimeStep > 0) {
                tmpMoleculeVelocityFixationList.add(
                    new MoleculeVelocityFixationDescription(
                        tmpMoleculeName, 
                        tmpIsFixedX, 
                        tmpIsFixedY, 
                        tmpIsFixedZ,
                        tmpVelocityX,
                        tmpVelocityY,
                        tmpVelocityZ,
                        tmpMaxTimeStep
                    )
                );
            }
        }
        if (!tmpMoleculeVelocityFixationList.isEmpty()) {
            return tmpMoleculeVelocityFixationList.toArray(new MoleculeVelocityFixationDescription[0]);
        } else {
            return null;
        }
    }
    
    /**
     * Molecule acceleration descriptions
     * 
     * @return Molecule acceleration descriptions or null if no molecule is accelerated
     */
    @Override
    public MoleculeAccelerationDescription[] getMoleculeAccelerationDescriptions() {
        String[][] tmpMoleculeAccelerationTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.MOLECULE_ACCELERATION);
        LinkedList<MoleculeAccelerationDescription> tmpMoleculeAccelerationList = new LinkedList<>();
        for (int i = 0; i < tmpMoleculeAccelerationTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeAccelerationTable[i][0];
            float tmpAccelerationX = Float.valueOf(tmpMoleculeAccelerationTable[i][1]);
            float tmpAccelerationY = Float.valueOf(tmpMoleculeAccelerationTable[i][2]);
            float tmpAccelerationZ = Float.valueOf(tmpMoleculeAccelerationTable[i][3]);
            int tmpFrequency = Integer.valueOf(tmpMoleculeAccelerationTable[i][4]);

            // If tmpMaxTimeStep is NOT defined it gets the highest possible value
            int tmpMaxTimeStep = Constants.MAXIMUM_NUMBER_OF_TIME_STEPS;
            if (tmpMoleculeAccelerationTable[i].length > 5) {
                tmpMaxTimeStep = Integer.valueOf(tmpMoleculeAccelerationTable[i][5]);
            }

            if ((tmpAccelerationX != 0.0f || tmpAccelerationY != 0.0f || tmpAccelerationZ != 0.0f) && tmpFrequency > 0  && tmpMaxTimeStep > 0) {
                tmpMoleculeAccelerationList.add(
                    new MoleculeAccelerationDescription(
                        tmpMoleculeName, 
                        tmpAccelerationX, 
                        tmpAccelerationY, 
                        tmpAccelerationZ,
                        tmpFrequency,
                        tmpMaxTimeStep
                    )
                );
            }
        }
        if (!tmpMoleculeAccelerationList.isEmpty()) {
            return tmpMoleculeAccelerationList.toArray(new MoleculeAccelerationDescription[0]);
        } else {
            return null;
        }
    }
    
    /**
     * Molecule particle radius of gyration (Rg) calculation descriptions
     * 
     * @return Molecule particle radius of gyration (Rg) calculation descriptions
     */
    @Override
    public MoleculeRgCalculationDescription[] getMoleculeRgCalculationDescriptions() {
        String[][] tmpMoleculeRgCalculationTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.RADIUS_OF_GYRATION);
        LinkedList<MoleculeRgCalculationDescription> tmpMoleculeRgCalculationList = new LinkedList<>();
        for (int i = 0; i < tmpMoleculeRgCalculationTable.length; i++) {
            String tmpMoleculeName = tmpMoleculeRgCalculationTable[i][0];
            boolean tmpIsRgCalculation = Boolean.valueOf(tmpMoleculeRgCalculationTable[i][1]);
            if (tmpIsRgCalculation) {
                tmpMoleculeRgCalculationList.add(new MoleculeRgCalculationDescription(tmpMoleculeName));
            }
        }
        if (!tmpMoleculeRgCalculationList.isEmpty()) {
            return tmpMoleculeRgCalculationList.toArray(new MoleculeRgCalculationDescription[0]);
        } else {
            return null;
        }
    }
    
    /**
     * Nearest-neighbor particle descriptions
     * 
     * @return Nearest-neighbor particle descriptions
     */
    @Override
    public NearestNeighborBaseParticleDescription[] getNearestNeighborParticleDescriptions() {
        String[][] tmpNearestNeighborParticleTable = this.getValueTable(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.NEAREST_NEIGHBOR_PARTICLE);
        LinkedList<NearestNeighborBaseParticleDescription> tmpNearestNeighborParticleList = new LinkedList<>();
        for (int i = 0; i < tmpNearestNeighborParticleTable.length; i++) {
            String tmpMoleculeName = tmpNearestNeighborParticleTable[i][0];
            String tmpParticleToken = tmpNearestNeighborParticleTable[i][1];
            boolean tmpIsNearestNeighborCalculation = Boolean.valueOf(tmpNearestNeighborParticleTable[i][2]);
            if (tmpIsNearestNeighborCalculation) {
                tmpNearestNeighborParticleList.add(new NearestNeighborBaseParticleDescription(tmpParticleToken, tmpMoleculeName));
            }
        }
        if (!tmpNearestNeighborParticleList.isEmpty()) {
            return tmpNearestNeighborParticleList.toArray(new NearestNeighborBaseParticleDescription[0]);
        } else {
            return null;
        }
    }
    
    /**
     * Nearest-neighbor distance in DPD units
     * 
     * @return Nearest-neighbor distance in DPD units
     */
    @Override
    public float getNearestNeighborDistance() {
        return this.getSingleFloatValue(Section.CHEMICAL_SYSTEM_DESCRIPTION, FileInput.NEAREST_NEIGHBOR_DISTANCE);
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Interaction description properties (get)">
    /**
     * Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * 
     * @return Temperature
     */
    @Override
    public float getTemperature() {
        return this.getSingleFloatValue(Section.INTERACTION_DESCRIPTION, FileInput.TEMPERATURE);
    }

    /**
     * DPD sigma parameter in DPD units
     * 
     * @return DPD sigma parameter
     */
    @Override
    public float getDpdSigma() {
        return this.getSingleFloatValue(Section.INTERACTION_DESCRIPTION, FileInput.DPD_SIGMA);
    }

    /**
     * Conservative force parameters a(ij)
     * 
     * @param aParticleTypes Particle types
     * @return Conservative force parameters a(ij)
     */
    @Override
    public float[][] getAij(ParticleTypes aParticleTypes) {
        String[][] tmpInteractionTable = this.getValueTable(Section.INTERACTION_DESCRIPTION, FileInput.INTERACTION_TABLE);
        // Quadratic matrix
        float[][] tmpAij = new float[aParticleTypes.getParticleTypeNumber()][aParticleTypes.getParticleTypeNumber()];
        for (String[] tmpInteractionTableLine : tmpInteractionTable) {
            int tmpIndexI = aParticleTypes.getIndex(tmpInteractionTableLine[0]);
            int tmpIndexJ = aParticleTypes.getIndex(tmpInteractionTableLine[1]);
            float tmpAijValue = Float.valueOf(tmpInteractionTableLine[2]);
            tmpAij[tmpIndexI][tmpIndexJ] = tmpAijValue;
            tmpAij[tmpIndexJ][tmpIndexI] = tmpAijValue;
        }
        return tmpAij;
    }
    
    /**
     * Hash map that maps particle-pair key to corresponding particle-pair bond.
     * Note: Particle-pair key has to be created with Utils.getParticlePairKey().
     * 
     * @return Hash map that maps particle-pair key to corresponding particle-pair bond or null if no bonds are defined
     */
    @Override
    public HashMap<String, ParticlePairHarmonicBond> getParticlePairBondMap() {
        if (this.hasParameter(Section.INTERACTION_DESCRIPTION, FileInput.BONDS_12_TABLE)) {
            String[][] tmpBonds12Table = this.getValueTable(Section.INTERACTION_DESCRIPTION, FileInput.BONDS_12_TABLE);
            HashMap<String, ParticlePairHarmonicBond> tmpParticlePairBondMap = new HashMap<>(tmpBonds12Table.length);
            for (String[] tmpBonds12Row : tmpBonds12Table) {
                String tmpParticlePairKey = Utils.getParticleTokenPairKey(tmpBonds12Row[0], tmpBonds12Row[1]);
                tmpParticlePairBondMap.put(
                    tmpParticlePairKey,
                    new ParticlePairHarmonicBond(
                        tmpParticlePairKey, 
                        Float.valueOf(tmpBonds12Row[2]), 
                        Float.valueOf(tmpBonds12Row[3]),
                        (Boolean.valueOf(tmpBonds12Row[4])) ? HarmonicBond.HarmonicBondBehaviour.DEFAULT : HarmonicBond.HarmonicBondBehaviour.ATTRACTIVE
                    )
                );
            }
            return tmpParticlePairBondMap;
        } else {
            return null;
        }
    }
    
    /**
     * Flag for use of Gaussian random for random force
     * 
     * @return True : Random DPD force is driven by random variable with 
     *                Gaussian distribution with zero mean and unit variance
     *                (slower)
     *         False: Random DPD force is driven by random variable with uniform
     *                distribution with zero mean and unit variance (faster)
     */
    @Override
    public boolean isGaussianRandomDpdForce() {
        return this.getSingleBooleanValue(Section.INTERACTION_DESCRIPTION, FileInput.IS_GAUSSIAN_RANDOM_DPD_FORCE);
    }
    
    /**
     * Electrostatics parameters
     * 
     * @return Electrostatics parameters or null if none are available
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    @Override
    public Electrostatics getElectrostatics() {
        if (this.hasParameter(Section.INTERACTION_DESCRIPTION, FileInput.ELECTROSTATICS)) {
            if (this.isElectrostaticsTypeDpd()) {
                // Factory.ElectrostaticsType.DPD
                String[] tmpElectrostaticsParameters = this.getNinetupelStringValues(Section.INTERACTION_DESCRIPTION, FileInput.ELECTROSTATICS);
                // Charge distribution
                Electrostatics.ChargeDistributionType tmpChargeDistributionType;
                if (tmpElectrostaticsParameters[6].equalsIgnoreCase(Electrostatics.ChargeDistributionType.NONE.toString())) {
                    tmpChargeDistributionType = Electrostatics.ChargeDistributionType.NONE;
                } else if (tmpElectrostaticsParameters[6].equalsIgnoreCase(Electrostatics.ChargeDistributionType.ALEJANDRE.toString())) {
                    tmpChargeDistributionType = Electrostatics.ChargeDistributionType.ALEJANDRE;
                } else {
                    throw new IllegalArgumentException("FileInput.getElectrostatics: Unknown charge distribution type.");
                }
                // Splitting type
                Electrostatics.SplittingType tmpSplittingType;
                if (tmpElectrostaticsParameters[8].equalsIgnoreCase(Electrostatics.SplittingType.NONE.toString())) {
                    tmpSplittingType = Electrostatics.SplittingType.NONE;
                } else if (tmpElectrostaticsParameters[8].equalsIgnoreCase(Electrostatics.SplittingType.FANOURGAKIS.toString())) {
                    tmpSplittingType = Electrostatics.SplittingType.FANOURGAKIS;
                } else {
                    throw new IllegalArgumentException("FileInput.getElectrostatics: Unknown splitting type.");
                }
                return
                    new Electrostatics(
                        Float.valueOf(tmpElectrostaticsParameters[0]), // aCutOffLength
                        Float.valueOf(tmpElectrostaticsParameters[1]), // aMaximumAbsoluteForceValue
                        Float.valueOf(tmpElectrostaticsParameters[2]), // anEffectiveExponent
                        Float.valueOf(tmpElectrostaticsParameters[3]), // aDampingDistance
                        Float.valueOf(tmpElectrostaticsParameters[4]), // aDampingFactor
                        Float.valueOf(tmpElectrostaticsParameters[5]), // anElectrostaticsCoupling
                        tmpChargeDistributionType,                      // aChargeDistributionType
                        Float.valueOf(tmpElectrostaticsParameters[7]), // aDecayLengthAlejandre
                        tmpSplittingType                                // aSplittingType
                    );
            } else {
                // Factory.ElectrostaticsType.AD_HOC
                float[] tmpElectrostaticsParameters = this.getSixtupelFloatValues(Section.INTERACTION_DESCRIPTION, FileInput.ELECTROSTATICS);
                return
                    new Electrostatics(
                        tmpElectrostaticsParameters[0], // aCutOffLength
                        tmpElectrostaticsParameters[1], // aMaximumAbsoluteForceValue
                        tmpElectrostaticsParameters[2], // anEffectiveChargeFactor
                        tmpElectrostaticsParameters[3], // anEffectiveExponent
                        tmpElectrostaticsParameters[4], // aDampingDistance
                        tmpElectrostaticsParameters[5]  // aDampingFactor
                    );
            }
        } else {
            return null;
        }
    }
    
    /**
     * Gravitational acceleration
     * 
     * @return Gravitational acceleration
     */
    @Override
    public GravitationalAcceleration getGravitationalAcceleration() {
        float[] tmpGraviationalAccelerationParameters =  this.getTripleFloatValues(Section.INTERACTION_DESCRIPTION, FileInput.GRAVITATIONAL_ACCELERATION);
        return new GravitationalAcceleration(tmpGraviationalAccelerationParameters[0], tmpGraviationalAccelerationParameters[1], tmpGraviationalAccelerationParameters[2]);
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Simulation description properties (get)">
    /**
     * Time step number for central simulation loop
     * 
     * @return Time step number
     */
    @Override
    public int getTimeStepNumber() {
        return this.getSingleIntegerValue(Section.SIMULATION_DESCRIPTION, FileInput.TIME_STEP_NUMBER);
    }

    /**
     * Time step length in DPD units
     * 
     * @return Time step length
     */
    @Override
    public float getTimeStepLength() {
        return this.getSingleFloatValue(Section.SIMULATION_DESCRIPTION, FileInput.TIME_STEP_LENGTH);
    }

    /**
     * Time step frequency for output
     * 
     * @return Time step frequency for output
     */
    @Override
    public int getTimeStepFrequencyForOutput() {
        return this.getSingleIntegerValue(Section.SIMULATION_DESCRIPTION, FileInput.TIME_STEP_FREQUENCY_FOR_OUTPUT);
    }
    
    /**
     * Number of initial potential energy minimization steps
     * 
     * @return Number of initial potential energy minimization steps
     */
    @Override
    public int getInitialPotentialEnergyMinimizationStepNumber() {
        return this.getFirstIntegerValue(Section.SIMULATION_DESCRIPTION, FileInput.INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_NUMBER);
    }
    
    /**
     * Type of initial potential energy minimization (true: All forces, false:
     * DPD force only)
     * 
     * @return Type of initial potential energy minimization
     */
    @Override
    public boolean isInitialPotentialEnergyMinimizationWithAllForces() {
        if (this.hasSecondBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_NUMBER)) {
            return this.getSecondBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_NUMBER);
        } else {
            // true: All forces for initial potential energy minimization
            return true;
        }
    }
    
    /**
     * Flag for initial potential energy minimization step output
     * 
     * @return True: Potential energy minimization step output is generated, false: Otherwise (NO output)
     */
    @Override
    public boolean isInitialPotentialEnergyMinimizationStepOutput() {
        return this.getSingleBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.IS_INITIAL_POTENTIAL_ENERGY_MINIMIZATION_STEP_OUTPUT);
    }
    
    /**
     * Periodic boundaries
     * 
     * @return Periodic boundaries
     */
    @Override
    public PeriodicBoundaries getPeriodicBoundaries() {
        boolean[] tmpPeriodicBoundariesParameters =  this.getTripleBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.PERIODIC_BOUNDARIES);
        return new PeriodicBoundaries(tmpPeriodicBoundariesParameters[0], tmpPeriodicBoundariesParameters[1], tmpPeriodicBoundariesParameters[2]);
    }

    /**
     * Flag for use of DPD unit masses:
     * True : DPD masses of all particles are set to 1
     * False: The DPD mass of the most lightweight particle (often water) is set to 1. 
     *        The masses of all other particles are set in relation to their 
     *        molar mass ratios to the most lightweight particle.
     * 
     * @return Flag value for use of DPD unit masses
     */
    @Override
    public boolean isDpdUnitMass() {
        return this.getSingleBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.IS_DPD_UNIT_MASS);
    }

    /**
     * Number of initial velocity scaling steps
     * 
     * @return Number of initial velocity scaling steps
     */
    @Override
    public int getInitialVelocityScalingSteps() {
        int tmpNumberOfInitialVelocityScalingSteps;
        if (this.hasParameter(Section.SIMULATION_DESCRIPTION, FileInput.INITIAL_VELOCITY_SCALING_STEPS)) {
            tmpNumberOfInitialVelocityScalingSteps = this.getSingleIntegerValue(Section.SIMULATION_DESCRIPTION, FileInput.INITIAL_VELOCITY_SCALING_STEPS);
        } else {
            // Code for compatibility
            boolean tmpIsVelocityScaling = this.getSingleBooleanValue(Section.SIMULATION_DESCRIPTION, FileInput.IS_VELOCITY_SCALING);
            if (tmpIsVelocityScaling) {
                tmpNumberOfInitialVelocityScalingSteps = Constants.MAXIMUM_NUMBER_OF_TIME_STEPS;
            } else {
                tmpNumberOfInitialVelocityScalingSteps = 0;
            }
        }
        return tmpNumberOfInitialVelocityScalingSteps;
    }
    
    /**
     * Seed for random number generation
     * 
     * @return Seed for random number generation
     */
    @Override
    public int getRandomSeed() {
        if (this.hasParameter(Section.SIMULATION_DESCRIPTION, FileInput.RANDOM_NUMBER_GENERATOR)) {
            String[] tmpValues = this.getTripleStringValues(Section.SIMULATION_DESCRIPTION, FileInput.RANDOM_NUMBER_GENERATOR);
            return Integer.valueOf(tmpValues[1]);
        } else {
            return FileInput.DEFAULT_RANDOM_SEED;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Simulation counts properties (get)">
    /**
     * Particle number
     * 
     * @return Particle number
     */
    @Override
    public int getParticleNumber() {
        return this.getSingleIntegerValue(Section.SIMULATION_COUNTS, FileInput.PARTICLE_NUMBER);
    }
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get/set)">
    /**
     * Maximum number of position correction trials
     * 
     * @return Maximum number of position correction trials
     */
    @Override
    public int getMaximumNumberOfPositionCorrectionTrials() {
        return this.maximumNumberOfPositionCorrectionTrials;
    }

    /**
     * Maximum number of position correction trials
     * 
     * @param aValue New value (only accepted if greater zero)
     */
    public void setMaximumNumberOfPositionCorrectionTrials(int aValue) {
        if (aValue > 0) {
            this.maximumNumberOfPositionCorrectionTrials = aValue;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Reads input file
     * 
     * @return True: Operation successful, false: Otherwise
     */
    private boolean readInputFile() {
        try {
            // <editor-fold defaultstate="collapsed" desc="Read sections">
            // Section GENERAL
            String[][] tmpSectionGeneral = this.readSectionFromInputFile(this.inputFilePathname, Section.GENERAL.toString());
            this.sectionMap.put(Section.GENERAL, tmpSectionGeneral);
            // Section PARTICLE_DESCRIPTION
            String[][] tmpSectionParticleDescription = this.readSectionFromInputFile(this.inputFilePathname, Section.PARTICLE_DESCRIPTION.toString());
            this.sectionMap.put(Section.PARTICLE_DESCRIPTION, tmpSectionParticleDescription);
            // Section CHEMICAL_SYSTEM_DESCRIPTION
            String[][] tmpSectionChemicalSystemDescription = this.readSectionFromInputFile(this.inputFilePathname, Section.CHEMICAL_SYSTEM_DESCRIPTION.toString());
            this.sectionMap.put(Section.CHEMICAL_SYSTEM_DESCRIPTION, tmpSectionChemicalSystemDescription);
            // Section INTERACTION_DESCRIPTION
            String[][] tmpSectionInteractionDescription = this.readSectionFromInputFile(this.inputFilePathname, Section.INTERACTION_DESCRIPTION.toString());
            this.sectionMap.put(Section.INTERACTION_DESCRIPTION, tmpSectionInteractionDescription);
            // Section SIMULATION_DESCRIPTION
            String[][] tmpSectionSimulationParameters = this.readSectionFromInputFile(this.inputFilePathname, Section.SIMULATION_DESCRIPTION.toString());
            this.sectionMap.put(Section.SIMULATION_DESCRIPTION, tmpSectionSimulationParameters);
            // Section SIMULATION_COUNTS
            String[][] tmpSectionCounts = this.readSectionFromInputFile(this.inputFilePathname, Section.SIMULATION_COUNTS.toString());
            this.sectionMap.put(Section.SIMULATION_COUNTS, tmpSectionCounts);
            // </editor-fold>
            return true;
        } catch (Exception anException) {
            return false;
        }
    }
    
    /**
     * Reads section from input file. Each line is split after one or more 
     * whitespace characters.
     *
     * @param anInputFilePathname Full pathname of input file (may be null then 
     * null is returned)
     * @param aSectionString Section string (if null/empty then null is returned)
     * @return Jagged string array from specified section of file or null if
     * jagged string array could not be read
     */
    private String[][] readSectionFromInputFile(String anInputFilePathname, String aSectionString) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anInputFilePathname == null || anInputFilePathname.isEmpty() || !(new File(anInputFilePathname)).isFile()) {
            return null;
        }
        if (aSectionString == null || aSectionString.isEmpty()) {
            return null;
        }
        // </editor-fold>
        BufferedReader tmpBufferedReader = null;
        try {
            FileReader tmpFileReader = new FileReader(anInputFilePathname);
            tmpBufferedReader = new BufferedReader(tmpFileReader, Constants.BUFFER_SIZE);
            LinkedList<String[]> tmpLinkedList = new LinkedList<>();
            String tmpLine;
            boolean tmpIsStarted = false;
            String tmpStartLine = String.format(FileInput.SECTION_START_TAG_FORMAT, aSectionString);
            String tmpEndLine = String.format(FileInput.SECTION_END_TAG_FORMAT, aSectionString);
            if (FileInput.COMMENT_LINE_PREFIX == null || FileInput.COMMENT_LINE_PREFIX.isEmpty()) {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    if (tmpIsStarted) {
                        if (tmpLine.trim().equalsIgnoreCase(tmpEndLine)) {
                            break;
                        }
                        String[] tmpItems = Utils.splitAndTrim(tmpLine.trim());
                        if (tmpItems != null) {
                            tmpLinkedList.add(tmpItems);
                        }
                    } else {
                        tmpIsStarted = tmpLine.trim().equalsIgnoreCase(tmpStartLine);
                    }
                }
            } else {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    if (tmpIsStarted) {
                        if (tmpLine.trim().equalsIgnoreCase(tmpEndLine)) {
                            break;
                        }
                        if (!tmpLine.startsWith(FileInput.COMMENT_LINE_PREFIX)) {
                            String[] tmpItems = Utils.splitAndTrim(tmpLine.trim());
                            if (tmpItems != null) {
                                tmpLinkedList.add(tmpItems);
                            }
                        }
                    } else {
                        tmpIsStarted = tmpLine.trim().equalsIgnoreCase(tmpStartLine);
                    }
                }
            }
            if (tmpLinkedList.size() > 0) {
                return tmpLinkedList.toArray(new String[0][]);
            } else {
                return null;
            }
        } catch (Exception e) {
            return null;
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException e) {
                    return null;
                }
            }
        }
    }

    /**
     * Reads positions and bonds file. Each line is split after one or more 
     * whitespace characters.
     *
     * @param aPositionsBondsFilePathname Full pathname of positions and bonds file (may be null then 
     * null is returned)
     * @return Jagged string array from file or null if jagged string array could not be read
     */
    private String[][] readPositionsBondsFile(String aPositionsBondsFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aPositionsBondsFilePathname == null || aPositionsBondsFilePathname.isEmpty() || !(new File(aPositionsBondsFilePathname)).isFile()) {
            return null;
        }
        // </editor-fold>
        BufferedReader tmpBufferedReader = null;
        try {
            FileReader tmpFileReader = new FileReader(aPositionsBondsFilePathname);
            tmpBufferedReader = new BufferedReader(tmpFileReader, Constants.BUFFER_SIZE);
            LinkedList<String[]> tmpLinkedList = new LinkedList<>();
            String tmpLine;
            if (FileInput.COMMENT_LINE_PREFIX == null || FileInput.COMMENT_LINE_PREFIX.isEmpty()) {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    String[] tmpItems = Utils.splitAndTrim(tmpLine.trim());
                    if (tmpItems != null) {
                        tmpLinkedList.add(tmpItems);
                    }
                }
            } else {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    if (!tmpLine.startsWith(FileInput.COMMENT_LINE_PREFIX)) {
                        String[] tmpItems = Utils.splitAndTrim(tmpLine.trim());
                        if (tmpItems != null) {
                            tmpLinkedList.add(tmpItems);
                        }
                    }
                }
            }
            if (tmpLinkedList.size() > 0) {
                return tmpLinkedList.toArray(new String[0][]);
            } else {
                return null;
            }
        } catch (Exception e) {
            return null;
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException e) {
                    return null;
                }
            }
        }
    }
    
    /**
     * Returns if section contains parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return True: Section contains parameter, false: Otherwise
     */
    private boolean hasParameter(Section aSection, String aParameterString) {
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    return true;
                }
            }
            return false;
        } catch (Exception anException) {
            return false;
        }
    }
    
    /**
     * Returns single string value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Single string value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private String getSingleStringValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSingleStringValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSingleStringValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 2) {
                        throw new IllegalArgumentException("FileInput.getSingleStringValue: tmpSectionLine has wrong format.");
                    }
                    return tmpSectionLine[1];
                }
            }
            throw new IllegalArgumentException("FileInput.getSingleStringValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSingleStringValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns multiple string values of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Multiple string values of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private String[] getMultipleStringValues(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getMultipleStringValues: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getMultipleStringValues: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length == 1) {
                        throw new IllegalArgumentException("FileInput.getMultipleStringValues: aParameterString does not have values.");
                    } else {
                        String[] tmpValues = new String[tmpSectionLine.length - 1];
                        for (int i = 1; i < tmpSectionLine.length; i++) {
                            tmpValues[i - 1] = tmpSectionLine[i];
                        }
                        return tmpValues;
                    }
                }
            }
            throw new IllegalArgumentException("FileInput.getMultipleStringValues: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getMultipleStringValues: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns single boolean value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Single boolean value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private boolean getSingleBooleanValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSingleBooleanValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSingleBooleanValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 2) {
                        throw new IllegalArgumentException("FileInput.getSingleBooleanValue: tmpSectionLine has wrong format.");
                    }
                    return Boolean.valueOf(tmpSectionLine[1]);
                }
            }
            throw new IllegalArgumentException("FileInput.getSingleBooleanValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSingleBooleanValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns second boolean value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Second boolean value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private boolean getSecondBooleanValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSecondBooleanValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSecondBooleanValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length < 3) {
                        throw new IllegalArgumentException("FileInput.getSecondBooleanValue: tmpSectionLine has wrong format.");
                    }
                    return Boolean.valueOf(tmpSectionLine[2]);
                }
            }
            throw new IllegalArgumentException("FileInput.getSecondBooleanValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSecondBooleanValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns if second boolean value of parameter exists
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return True: Second boolean value of parameter exists, false: Otherwise
     */
    private boolean hasSecondBooleanValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            return false;
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            return false;
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length < 3) {
                        return false;
                    }
                    boolean tmpSecondBooleanValue = Boolean.valueOf(tmpSectionLine[2]);
                    return true;
                }
            }
            return false;
        } catch (Exception anException) {
            return false;
        }
    }
    
    /**
     * Returns triple boolean value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Triple boolean value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private boolean[] getTripleBooleanValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getTripleBooleanValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getTripleBooleanValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 4) {
                        throw new IllegalArgumentException("FileInput.getTripleBooleanValue: tmpSectionLine has wrong format.");
                    }
                    return new boolean[] {Boolean.valueOf(tmpSectionLine[1]), Boolean.valueOf(tmpSectionLine[2]), Boolean.valueOf(tmpSectionLine[3])};
                }
            }
            throw new IllegalArgumentException("FileInput.getTripleBooleanValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getTripleBooleanValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns single integer value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Single integer value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private int getSingleIntegerValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSingleIntegerValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSingleIntegerValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 2) {
                        throw new IllegalArgumentException("FileInput.getSingleIntegerValue: tmpSectionLine has wrong format.");
                    }
                    return Integer.valueOf(tmpSectionLine[1]);
                }
            }
            throw new IllegalArgumentException("FileInput.getSingleIntegerValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSingleIntegerValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns first integer value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return First integer value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private int getFirstIntegerValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getFirstIntegerValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getFirstIntegerValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length < 2) {
                        throw new IllegalArgumentException("FileInput.getFirstIntegerValue: tmpSectionLine has wrong format.");
                    }
                    return Integer.valueOf(tmpSectionLine[1]);
                }
            }
            throw new IllegalArgumentException("FileInput.getFirstIntegerValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getFirstIntegerValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns single float value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Single float value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private float getSingleFloatValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSingleFloatValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSingleFloatValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 2) {
                        throw new IllegalArgumentException("FileInput.getSingleFloatValue: tmpSectionLine has wrong format.");
                    }
                    return Float.valueOf(tmpSectionLine[1]);
                }
            }
            throw new IllegalArgumentException("FileInput.getSingleFloatValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSingleFloatValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns triple float value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Triple float value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private float[] getTripleFloatValues(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getTripleFloatValues: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getTripleFloatValues: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 4) {
                        throw new IllegalArgumentException("FileInput.getTripleFloatValues: tmpSectionLine has wrong format.");
                    }
                    return new float[] {Float.valueOf(tmpSectionLine[1]), Float.valueOf(tmpSectionLine[2]), Float.valueOf(tmpSectionLine[3])};
                }
            }
            throw new IllegalArgumentException("FileInput.getTripleFloatValues: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getTripleFloatValues: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns triple String value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Triple String value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private String[] getTripleStringValues(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getTripleStringValues: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getTripleStringValues: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 4) {
                        throw new IllegalArgumentException("FileInput.getTripleStringValues: tmpSectionLine has wrong format.");
                    }
                    return new String[] {tmpSectionLine[1], tmpSectionLine[2], tmpSectionLine[3]};
                }
            }
            throw new IllegalArgumentException("FileInput.getTripleStringValues: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getTripleStringValues: aParameterString can not be evaluated.", anException);
        }
    }

    /**
     * Returns ninetupel String value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Ninetupel String value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private String[] getNinetupelStringValues(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getNinetupelStringValues: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getNinetupelStringValues: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 10) {
                        throw new IllegalArgumentException("FileInput.getNinetupelStringValues: tmpSectionLine has wrong format.");
                    }
                    return new String[] {
                        tmpSectionLine[1], 
                        tmpSectionLine[2], 
                        tmpSectionLine[3],
                        tmpSectionLine[4],
                        tmpSectionLine[5],
                        tmpSectionLine[6],
                        tmpSectionLine[7],
                        tmpSectionLine[8],
                        tmpSectionLine[9]
                    };
                }
            }
            throw new IllegalArgumentException("FileInput.getNinetupelStringValues: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getNinetupelStringValues: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns sixtupel float value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Sixtupel float  value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private float[] getSixtupelFloatValues(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSixtupelFloatValues: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSixtupelFloatValues: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 7) {
                        throw new IllegalArgumentException("FileInput.getSixtupelFloatValues: tmpSectionLine has wrong format.");
                    }
                    return new float[] {
                        Float.valueOf(tmpSectionLine[1]), 
                        Float.valueOf(tmpSectionLine[2]), 
                        Float.valueOf(tmpSectionLine[3]),
                        Float.valueOf(tmpSectionLine[4]), 
                        Float.valueOf(tmpSectionLine[5]), 
                        Float.valueOf(tmpSectionLine[6])
                    };
                }
            }
            throw new IllegalArgumentException("FileInput.getSixtupelFloatValues: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSixtupelFloatValues: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns single long value of parameter
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Single long value of parameter
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private long getSingleLongValue(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getSingleLongValue: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getSingleLongValue: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (String[] tmpSectionLine : tmpSection) {
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (tmpSectionLine.length != 2) {
                        throw new IllegalArgumentException("FileInput.getSingleLongValue: tmpSectionLine has wrong format.");
                    }
                    return Long.valueOf(tmpSectionLine[1]);
                }
            }
            throw new IllegalArgumentException("FileInput.getSingleLongValue: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getSingleLongValue: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns value table
     * 
     * @param aSection Section tag
     * @param aParameterString Parameter string
     * @return Value table
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    private String[][] getValueTable(Section aSection, String aParameterString) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSection == null || !this.sectionMap.containsKey(aSection)) {
            throw new IllegalArgumentException("FileInput.getParameterTable: aSection is null/empty or does not exist.");
        }
        if (aParameterString == null || aParameterString.isEmpty()) {
            throw new IllegalArgumentException("FileInput.getParameterTable: aParameterString is null/empty.");
        }
        // </editor-fold>
        try {
            String[][] tmpSection = this.sectionMap.get(aSection);
            for (int i = 0; i < tmpSection.length; i++) {
                String[] tmpSectionLine = tmpSection[i];
                if (tmpSectionLine[0].equalsIgnoreCase(aParameterString)) {
                    if (!tmpSection[i + 1][0].equals(FileInput.TABLE_START)) {
                        throw new IllegalArgumentException("FileInput.getParameterTable: No TABLE_START found.");
                    }
                    LinkedList<String[]> tmpLineList = new LinkedList<>();
                    for (int k = i + 2; k < tmpSection.length; k++) {
                        tmpSectionLine = tmpSection[k];
                        if (tmpSectionLine[0].equals(FileInput.TABLE_END)) {
                            return tmpLineList.toArray(new String[0][]);
                        } else {
                            tmpLineList.add(tmpSectionLine);
                        }
                    }
                    throw new IllegalArgumentException("FileInput.getParameterTable: No TABLE_END found.");
                }
            }
            throw new IllegalArgumentException("FileInput.getParameterTable: aParameterString could not be found in section.");
        } catch (Exception anException) {
            throw new IllegalArgumentException("FileInput.getParameterTable: aParameterString can not be evaluated.", anException);
        }
    }
    
    /**
     * Returns string reference to particle token
     * 
     * @param aParticleToken Particle token
     * @return String reference to particle token
     */
    private String getReferencedParticleToken(String aParticleToken) {
        String tmpReferencedParticleToken = this.particleTokenMap.get(aParticleToken);
        if (tmpReferencedParticleToken == null) {
            this.particleTokenMap.put(aParticleToken, aParticleToken);
            tmpReferencedParticleToken = aParticleToken;
        }
        return tmpReferencedParticleToken;
    }
    
    /**
     * Returns if electrostatics type is DPD
     * 
     * @return True: Electrostatics type is DPD, false: Otherwise
     */
    private boolean isElectrostaticsTypeDpd() {
        String[][] tmpSection = this.sectionMap.get(Section.INTERACTION_DESCRIPTION);
        for (String[] tmpSectionLine : tmpSection) {
            if (tmpSectionLine[0].equalsIgnoreCase(FileInput.ELECTROSTATICS)) {
                if (tmpSectionLine.length != 7) {
                    // Factory.ElectrostaticsType.DPD
                    return true;
                } else {
                    return false;
                }
            }
        }
        return false;
    }
    // </editor-fold>
    
}
