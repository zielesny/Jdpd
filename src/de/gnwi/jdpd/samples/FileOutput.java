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
package de.gnwi.jdpd.samples;

import de.gnwi.jdpd.interfaces.IOutput;
import de.gnwi.jdpd.parameters.RestartInfo;
import de.gnwi.jdpd.particlePosition.ParticlePosition;
import de.gnwi.jdpd.particlePosition.ParticlePositionPool;
import de.gnwi.jdpd.rg.MoleculeRgValue;
import de.gnwi.jdpd.utilities.Constants;
import de.gnwi.jdpd.utilities.FileOutputStrings;
import de.gnwi.jdpd.utilities.Strings;
import de.gnwi.jdpd.utilities.Utils;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.GZIPOutputStream;

/**
 * File output
 * 
 * @author Achim Zielesny
 */
public class FileOutput implements IOutput {

    // <editor-fold defaultstate="collapsed" desc="Private class ParticlePositionsWriteTask">
    /**
     * Runnable for writing particle positions to file in a parallelised manner with 
     * executor service.
     */
    private class ParticlePositionsWriteTask implements Runnable {
        
        // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
        /**
         * Version
         */
        private final static String VERSION_1_0_0 = "Version 1.0.0";
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Private final class variables">
        /**
         * Particle positions file pathname
         */
        private final String particlePositionsFilePathname;
        
        /**
         * Particle positions
         */
        private final ParticlePosition[] particlePositions;
        
        /**
         * Queue with particle positions file pathnames
         */
        private final ConcurrentLinkedQueue<String> particlePositionsFilePathnameQueue;
    
        /**
         * Particle position pools
         */
        private final ParticlePositionPool particlePositionPool;
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Private class variables">
        /**
         * Number of after-decimal-separator digits for particle positions (unrestricted if smaller of equal to 0)
         */
        private final int numberOfAfterDecimalDigitsForParticlePositions;
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Constructor">
        /**
         * Constructor
         * Note: NO checks are performed.
         * 
         * @param aParticlePositionsFilePathname Particle positions file pathname
         * @param aParticlePositions Particle positions 
         * @param aParticlePositionsFilePathnameQueue Queue with particle positions step file pathnames (may be null)
         * @param aParticlePositionPool Particle position pool (may be null)
         * @param aNumberOfAfterDecimalDigitsForParticlePositions Number of after-decimal-separator digits for particle positions (unrestricted if smaller of equal to 0)
         */
        public ParticlePositionsWriteTask(
            String aParticlePositionsFilePathname,
            ParticlePosition[] aParticlePositions,
            ConcurrentLinkedQueue<String> aParticlePositionsFilePathnameQueue,
            ParticlePositionPool aParticlePositionPool,
            int aNumberOfAfterDecimalDigitsForParticlePositions
            ) {
            this.particlePositionsFilePathname = aParticlePositionsFilePathname;
            this.particlePositions = aParticlePositions;
            this.particlePositionsFilePathnameQueue = aParticlePositionsFilePathnameQueue;
            this.particlePositionPool = aParticlePositionPool;
            this.numberOfAfterDecimalDigitsForParticlePositions = aNumberOfAfterDecimalDigitsForParticlePositions;
        }
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Public overridden methods">
        /**
         * Overridden run method
         * Note: Particle position instances are set in ParticlePositionPool 
         * for re-use.
         */
        @Override
        public void run() {
            // IMPORTANT: Sort particle positions first
            Arrays.sort(this.particlePositions);
            try (PrintWriter tmpPrintWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(this.particlePositionsFilePathname), Constants.BUFFER_SIZE));) {
                this.writeOutput(tmpPrintWriter);
            } catch (Exception e) {
                return;
            }
            // Set particle positions for re-use in ParticlePositionPool 
            if (this.particlePositionPool != null) {
                for (ParticlePosition tmpParticlePosition : this.particlePositions) {
                    this.particlePositionPool.setParticlePositionForReuse(tmpParticlePosition);
                }
            }
            // Add this.particlePositionsFilePathname to this.particlePositionsFilePathnameQueue if available
            this.particlePositionsFilePathnameQueue.add(this.particlePositionsFilePathname);
        }
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Private methods">
        /**
         * Writes output (see code)
         * NOTE: No checks are performed!
         * 
         * @param aPrintWriter Print writer to write with
         */
        private void writeOutput(PrintWriter aPrintWriter) {
            // <editor-fold defaultstate="collapsed" desc="Write version">
            aPrintWriter.println(VERSION_1_0_0);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Write number of particles">
            aPrintWriter.println(String.valueOf(this.particlePositions.length));
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Write particle positions">
            DecimalFormat tmpDecimalFormatter = this.getDecimalFormatter(this.numberOfAfterDecimalDigitsForParticlePositions);
            int tmpStartIndex = 0;
            for (int i = 0; i < this.particlePositions.length; i++) {
                String tmpComparativeMoleculeParticleString = null;
                int tmpCounter = 0;
                for (int k = tmpStartIndex; k < this.particlePositions.length; k++) {
                    if (tmpComparativeMoleculeParticleString == null) {
                        tmpComparativeMoleculeParticleString = this.particlePositions[k].getMoleculeParticleString();
                        tmpCounter++;
                    } else {
                        if (tmpComparativeMoleculeParticleString.equals(this.particlePositions[k].getMoleculeParticleString())) {
                            tmpCounter++;
                        } else {
                            break;
                        }
                    }
                }
                ParticlePosition tmpParticlePosition = this.particlePositions[tmpStartIndex];
                // Write molecule name
                aPrintWriter.println(tmpParticlePosition.getMoleculeName());
                // Write particle token
                aPrintWriter.println(tmpParticlePosition.getParticleToken());
                // Write number of x,y,z-positions
                aPrintWriter.println(String.valueOf(tmpCounter));
                // Write x,y,z-positions plus particle index and molecule index
                if (tmpDecimalFormatter == null) {
                    for (int k = 0; k < tmpCounter; k++) {
                        tmpParticlePosition = this.particlePositions[tmpStartIndex + k];
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getXPosition()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getYPosition()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getZPosition()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getParticleIndex()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getMoleculeIndex()));
                    }
                } else {
                    for (int k = 0; k < tmpCounter; k++) {
                        tmpParticlePosition = this.particlePositions[tmpStartIndex + k];
                        aPrintWriter.println(tmpDecimalFormatter.format(tmpParticlePosition.getXPosition()));
                        aPrintWriter.println(tmpDecimalFormatter.format(tmpParticlePosition.getYPosition()));
                        aPrintWriter.println(tmpDecimalFormatter.format(tmpParticlePosition.getZPosition()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getParticleIndex()));
                        aPrintWriter.println(String.valueOf(tmpParticlePosition.getMoleculeIndex()));
                    }
                }
                tmpStartIndex += tmpCounter;
            }
            // </editor-fold>
        }
        
        /**
         * Returns decimal formatter with specified number of 
         * after-decimal-separator digits.
         * 
         * @param aNumberOfAfterDecimalSeparatorDigits Number of 
         * after-decimal-separator digits
         * @return Decimal formatter with specified number of 
         * after-decimal-separator digits of null if number was smaller or equal 
         * to zero.
         */
        private DecimalFormat getDecimalFormatter(int aNumberOfAfterDecimalSeparatorDigits) {
            if (aNumberOfAfterDecimalSeparatorDigits > 0) {
                String tmpDecimalFormatString = "#." + String.join("", Collections.nCopies(aNumberOfAfterDecimalSeparatorDigits, "#"));
                DecimalFormat tmpDecimalFormat = new DecimalFormat(tmpDecimalFormatString);
                DecimalFormatSymbols tmpDecimalFormatSymbols = DecimalFormatSymbols.getInstance();
                // IMPORTANT: Set dot as decimal separator
                tmpDecimalFormatSymbols.setDecimalSeparator('.');
                tmpDecimalFormat.setDecimalFormatSymbols(tmpDecimalFormatSymbols);
                return tmpDecimalFormat;
            } else {
                return null;
            }
        }
        // </editor-fold>

    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Executor service
     */
    private final ExecutorService executorService;

    /**
     * Simulation step file pathname
     */
    private final String simulationStepFilePathname;

    /**
     * Temperature file pathname
     */
    private final String temperatureFilePathname;

    /**
     * DPD potential energy file pathname
     */
    private final String uPotDpdFilePathname;

    /**
     * Bond potential energy file pathname
     */
    private final String uPotBondFilePathname;

    /**
     * Electrostatics potential energy file pathname
     */
    private final String uPotElectrostaticsFilePathname;

    /**
     * Total potential energy file pathname
     */
    private final String uPotTotalFilePathname;

    /**
     * Kinetic energy file pathname
     */
    private final String uKinFilePathname;

    /**
     * Total energy file pathname
     */
    private final String uTotalFilePathname;

    /**
     * Surface tension along x axis file pathname
     */
    private final String surfaceTensionAlongXFilePathname;

    /**
     * Surface tension along y axis file pathname
     */
    private final String surfaceTensionAlongYFilePathname;

    /**
     * Surface tension along z axis file pathname
     */
    private final String surfaceTensionAlongZFilePathname;

    /**
     * Surface tension norm file pathname
     */
    private final String surfaceTensionNormFilePathname;

    /**
     * DPD surface tension along x axis file pathname
     */
    private final String dpdSurfaceTensionAlongXFilePathname;

    /**
     * DPD surface tension along y axis file pathname
     */
    private final String dpdSurfaceTensionAlongYFilePathname;

    /**
     * DPD surface tension along z axis file pathname
     */
    private final String dpdSurfaceTensionAlongZFilePathname;

    /**
     * DPD surface tension norm file pathname
     */
    private final String dpdSurfaceTensionNormFilePathname;
    
    /**
     * Output directory path
     */
    private final String outputDirectoryPath;

    /**
     * Radius-of-gyration (Rg) directory path
     */
    private final String radiusOfGyrationDirectoryPath;

    /**
     * File pathname for base molecule-particle to nearest-neighbor molecule-particle simulation step vs. frequency map
     */
    private final String baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname;

    /**
     * File pathname for base molecule-particle to nearest-neighbor particle simulation step vs. frequency map
     */
    private final String baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname;

    /**
     * File pathname for base molecule-particle to nearest-neighbor molecule simulation step vs. frequency map
     */
    private final String baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname;

    /**
     * File pathname for base molecule to nearest-neighbor molecule simulation step vs. frequency map
     */
    private final String baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname;

    /**
     * File pathname for base molecule to nearest-neighbor molecule-tuple simulation step vs. frequency map
     */
    private final String baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname;
    
    /**
     * Directory path for simulation step particle positions
     */
    private final String simulationStepParticlePositionsDirectoryPath;
    
    /**
     * Directory path for minimisation step particle positions
     */
    private final String minimizationStepParticlePositionsDirectoryPath;
    
    /**
     * Number of parallel tasks for write operations
     */
    private final int parallelTaskNumber;
    
    /**
     * Queue with particle positions file pathnames
     */
    private final ConcurrentLinkedQueue<String> particlePositionsfilePathnameQueue;
    
    /**
     * Simulation step list
     */
    private final LinkedList<String> simulationStepList;
    
    /**
     * Temperature list
     */
    private final LinkedList<String> temperatureList;
    
    /**
     * DPD potential energy list
     */
    private final LinkedList<String> uPotDpdList;
    
    /**
     * Bond potential energy list
     */
    private final LinkedList<String> uPotBondList;
    
    /**
     * Electrostatics potential energy list
     */
    private final LinkedList<String> uPotElectrostaticsList;
    
    /**
     * Total potential energy list
     */
    private final LinkedList<String> uPotTotalList;
    
    /**
     * Kinetic energy list
     */
    private final LinkedList<String> uKinList;
    
    /**
     * Total energy list
     */
    private final LinkedList<String> uTotalList;
    
    /**
     * Surface tension along x axis list
     */
    private final LinkedList<String> surfaceTensionAlongXList;
    
    /**
     * Surface tension along y axis list
     */
    private final LinkedList<String> surfaceTensionAlongYList;
    
    /**
     * Surface tension along z axis list
     */
    private final LinkedList<String> surfaceTensionAlongZList;
    
    /**
     * Surface tension norm list
     */
    private final LinkedList<String> surfaceTensionNormList;
    
    /**
     * DPD surface tension along x axis list
     */
    private final LinkedList<String> dpdSurfaceTensionAlongXList;
    
    /**
     * DPD surface tension along y axis list
     */
    private final LinkedList<String> dpdSurfaceTensionAlongYList;
    
    /**
     * DPD surface tension along z axis list
     */
    private final LinkedList<String> dpdSurfaceTensionAlongZList;
    
    /**
     * DPD surface tension norm list
     */
    private final LinkedList<String> dpdSurfaceTensionNormList;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Counter for setSimulationStepInformation() calls
     */
    private int setCount;
    
    /**
     * Last simulation step
     */
    private int lastSimulationStep;
    
    /**
     * File pathname of particle positions of last simulation step
     */
    private String lastSimulationStepParticlePositionsFilePathname;
    
    /**
     * Restart info
     */
    private RestartInfo restartInfo;
    
    /**
     * Number of after-decimal-separator digits for particle positions (unrestricted if smaller of equal to 0)
     */
    private int numberOfAfterDecimalDigitsForParticlePositions;
    
    /**
     * Particle position pools
     */
    private ParticlePositionPool particlePositionPool;
    
    /**
     * Molecule names of Rg calculation
     * NOTE: rgValueLists[i] corresponds to rgMoleculeNames[i]
     */
    private String[] rgMoleculeNames;

    /**
     * Lists with Rg values
     * NOTE: rgValueLists[i] corresponds to rgMoleculeNames[i]
     */
    private LinkedList<String>[] rgValueLists;

    /**
     * Base molecule-particle to nearest-neighbor molecule-particle simulation step vs. frequency map
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap;

    /**
     * Base molecule-particle to nearest-neighbor particle simulation step vs. frequency map
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap;

    /**
     * Base molecule-particle to nearest-neighbor molecule simulation step vs. frequency map
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap;

    /**
     * Base molecule to nearest-neighbor molecule simulation step vs. frequency map
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> baseMoleculeToNearestNeighborMoleculeStepFrequencyMap;

    /**
     * Base molecule to nearest-neighbor molecule-tuple simulation step vs. frequency map
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param anOutputDirectoryPath Output directory path
     * @param aPropertiesDirectoryPath Property directory path
     * @param aRadiusOfGyrationDirectoryPath Radius of gyration (Rg) directory path
     * @param aNearestNeighborDirectoryPath NearestNeighbor directory path
     * @param aSimulationStepParticlePositionsDirectoryPath Directory path for simulation step particle positions (may be null then nothing is written)
     * @param aMinimizationStepParticlePositionsDirectoryPath Directory path for minimisation step particle positions (may be null then nothing is written)
     * @param aParallelTaskNumber Number of parallel tasks for write operations
     * @param aNumberOfAfterDecimalDigitsForParticlePositions Number of after-decimal-separator digits for particle positions (unrestricted if smaller of equal to 0)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public FileOutput (
        String anOutputDirectoryPath,
        String aPropertiesDirectoryPath,
        String aRadiusOfGyrationDirectoryPath,
        String aNearestNeighborDirectoryPath,
        String aSimulationStepParticlePositionsDirectoryPath,
        String aMinimizationStepParticlePositionsDirectoryPath,
        int aParallelTaskNumber,
        int aNumberOfAfterDecimalDigitsForParticlePositions
    ) throws IllegalArgumentException 
    {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anOutputDirectoryPath == null || anOutputDirectoryPath.isEmpty()) {
            throw new IllegalArgumentException("FileOutput.Constructor: anOutputDirectoryPath is null/empty.");
        }
        if (aPropertiesDirectoryPath == null || aPropertiesDirectoryPath.isEmpty()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aPropertyDirectoryPath is null/empty.");
        }
        if (aParallelTaskNumber < 1) {
            throw new IllegalArgumentException("FileOutput.Constructor: aParallelTaskNumber is less than 1.");
        }
        if (!new File(anOutputDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: anOutputDirectoryPath does not exist.");
        }
        if (!new File(aPropertiesDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aPropertiesDirectoryPath does not exist.");
        }
        if (aRadiusOfGyrationDirectoryPath != null && !new File(aRadiusOfGyrationDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aRadiusOfGyrationDirectoryPath does not exist.");
        }
        if (aNearestNeighborDirectoryPath != null && !new File(aNearestNeighborDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aNearestNeighborDirectoryPath does not exist.");
        }
        if (aSimulationStepParticlePositionsDirectoryPath != null && !new File(aSimulationStepParticlePositionsDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aSimulationStepParticlePositionsDirectoryPath does not exist.");
        }
        if (aMinimizationStepParticlePositionsDirectoryPath != null && !new File(aMinimizationStepParticlePositionsDirectoryPath).isDirectory()) {
            throw new IllegalArgumentException("FileOutput.Constructor: aMinimizationStepParticlePositionsDirectoryPath does not exist.");
        }
        // </editor-fold>
        this.outputDirectoryPath = anOutputDirectoryPath;
        this.simulationStepFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.SIMULATION_STEP_FILENAME;
        this.temperatureFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.TEMPERATURE_FILENAME;
        this.uPotDpdFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_POT_DPD_FILENAME;
        this.uPotBondFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_POT_BOND_FILENAME;
        this.uPotElectrostaticsFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_POT_ELECTROSTATICS_FILENAME;
        this.uPotTotalFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_POT_TOTAL_FILENAME;
        this.uKinFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_KIN_FILENAME;
        this.uTotalFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.U_TOTAL_FILENAME;
        this.surfaceTensionAlongXFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.SURFACE_TENSION_X_FILENAME;
        this.surfaceTensionAlongYFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.SURFACE_TENSION_Y_FILENAME;
        this.surfaceTensionAlongZFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.SURFACE_TENSION_Z_FILENAME;
        this.surfaceTensionNormFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.SURFACE_TENSION_NORM_FILENAME;
        this.dpdSurfaceTensionAlongXFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.DPD_SURFACE_TENSION_X_FILENAME;
        this.dpdSurfaceTensionAlongYFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.DPD_SURFACE_TENSION_Y_FILENAME;
        this.dpdSurfaceTensionAlongZFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.DPD_SURFACE_TENSION_Z_FILENAME;
        this.dpdSurfaceTensionNormFilePathname = aPropertiesDirectoryPath + File.separatorChar + FileOutputStrings.DPD_SURFACE_TENSION_NORM_FILENAME;

        this.radiusOfGyrationDirectoryPath = aRadiusOfGyrationDirectoryPath;

        if (aNearestNeighborDirectoryPath != null) {
            this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname = 
                aNearestNeighborDirectoryPath + File.separatorChar + FileOutputStrings.MP_TO_NN_MP_FILENAME;
            this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname =
                aNearestNeighborDirectoryPath + File.separatorChar + FileOutputStrings.MP_TO_NN_P_FILENAME;
            this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname =
                aNearestNeighborDirectoryPath + File.separatorChar + FileOutputStrings.MP_TO_NN_M_FILENAME;
            this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname =
                aNearestNeighborDirectoryPath + File.separatorChar + FileOutputStrings.M_TO_NN_M_FILENAME;
            this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname =
                aNearestNeighborDirectoryPath + File.separatorChar + FileOutputStrings.M_TO_NN_M_TUPLE_FILENAME;
        } else {
            this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname = null;
            this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname = null;
            this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname = null;
            this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname = null;
            this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname = null;
        }
        
        this.simulationStepParticlePositionsDirectoryPath = aSimulationStepParticlePositionsDirectoryPath;
        this.minimizationStepParticlePositionsDirectoryPath = aMinimizationStepParticlePositionsDirectoryPath;
        
        this.parallelTaskNumber = aParallelTaskNumber;
        this.executorService = Executors.newFixedThreadPool(this.parallelTaskNumber);
        this.particlePositionsfilePathnameQueue = new ConcurrentLinkedQueue<>();
     
        this.numberOfAfterDecimalDigitsForParticlePositions = aNumberOfAfterDecimalDigitsForParticlePositions;
        this.particlePositionPool = null;
        
        this.simulationStepList = this.getInitializedList(this.simulationStepFilePathname);
        this.temperatureList = this.getInitializedList(this.temperatureFilePathname);
        this.uPotDpdList = this.getInitializedList(this.uPotDpdFilePathname);
        this.uPotBondList = this.getInitializedList(this.uPotBondFilePathname);
        this.uPotElectrostaticsList = this.getInitializedList(this.uPotElectrostaticsFilePathname);
        this.uPotTotalList = this.getInitializedList(this.uPotTotalFilePathname);
        this.uKinList = this.getInitializedList(this.uKinFilePathname);
        this.uTotalList = this.getInitializedList(this.uTotalFilePathname);
        this.surfaceTensionAlongXList = this.getInitializedList(this.surfaceTensionAlongXFilePathname);
        this.surfaceTensionAlongYList = this.getInitializedList(this.surfaceTensionAlongYFilePathname);
        this.surfaceTensionAlongZList = this.getInitializedList(this.surfaceTensionAlongZFilePathname);
        this.surfaceTensionNormList = this.getInitializedList(this.surfaceTensionNormFilePathname);
        this.dpdSurfaceTensionAlongXList = this.getInitializedList(this.dpdSurfaceTensionAlongXFilePathname);
        this.dpdSurfaceTensionAlongYList = this.getInitializedList(this.dpdSurfaceTensionAlongYFilePathname);
        this.dpdSurfaceTensionAlongZList = this.getInitializedList(this.dpdSurfaceTensionAlongZFilePathname);
        this.dpdSurfaceTensionNormList = this.getInitializedList(this.dpdSurfaceTensionNormFilePathname);
        
        this.rgMoleculeNames = null;
        this.rgValueLists = null;
        
        this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap = null;
        this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap = null;
        this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap = null;
        this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMap = null;
        this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap = null;
        
        this.setCount = 0;
        this.lastSimulationStep = -1;
        this.lastSimulationStepParticlePositionsFilePathname = null;
        
        this.restartInfo = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Sets particle positions at start
     * 
     * @param aParticlePositions Particle positions
     */
    @Override
    public void setStartParticlePositions(ParticlePosition[] aParticlePositions) {
        this.setCount++;
        String tmpFileEnding;
        tmpFileEnding = Strings.GZIP_FILE_ENDING;
        String tmpParticlePositionsStartFilePathname = this.outputDirectoryPath + File.separatorChar + FileOutputStrings.PARTICLE_POSITIONS_START_FILENAME_PREFIX + tmpFileEnding;
        this.executorService.execute(
            new ParticlePositionsWriteTask(
                tmpParticlePositionsStartFilePathname, 
                aParticlePositions, 
                this.particlePositionsfilePathnameQueue,
                this.particlePositionPool,
                this.numberOfAfterDecimalDigitsForParticlePositions
            )
        );
    }
    
    /**
     * Sets particle positions of minimisation step
     * 
     * @param aMinimizationStep Minimisation step
     * @param aParticlePositions Particle positions
     */
    @Override
    public void setMinimizationStepParticlePositions(int aMinimizationStep, ParticlePosition[] aParticlePositions) {
        if (this.minimizationStepParticlePositionsDirectoryPath != null) {
            this.setCount++;
            String tmpFileEnding;
            tmpFileEnding = Strings.GZIP_FILE_ENDING;
            String tmpParticlePositionsMinimizationStepFilePathname = String.format(FileOutputStrings.PARTICLE_POSITIONS_MINIMIZATION_STEP_FILE_PREFIX_FORMAT, this.minimizationStepParticlePositionsDirectoryPath + File.separatorChar, String.valueOf(aMinimizationStep)) + tmpFileEnding;
            this.executorService.execute(
                new ParticlePositionsWriteTask(
                    tmpParticlePositionsMinimizationStepFilePathname, 
                    aParticlePositions, 
                    this.particlePositionsfilePathnameQueue, 
                    this.particlePositionPool,
                    this.numberOfAfterDecimalDigitsForParticlePositions
                )
            );
        }
    }
    
    /**
     * Sets particle positions after optimisation
     * 
     * @param aParticlePositions Particle positions
     */
    @Override
    public void setMinimizedParticlePositions(ParticlePosition[] aParticlePositions) {
        this.setCount++;
        String tmpFileEnding;
        tmpFileEnding = Strings.GZIP_FILE_ENDING;
        String tmpParticlePositionsMinimizedFilePathname = this.outputDirectoryPath + File.separatorChar + FileOutputStrings.PARTICLE_POSITIONS_MINIMIZED_FILENAME_PREFIX + tmpFileEnding;
        this.executorService.execute(
            new ParticlePositionsWriteTask(
                tmpParticlePositionsMinimizedFilePathname, 
                aParticlePositions, 
                this.particlePositionsfilePathnameQueue, 
                this.particlePositionPool,
                this.numberOfAfterDecimalDigitsForParticlePositions
            )
        );
    }

    /**
     * Sets simulation step information
     * Note: NO checks are performed.
     * Note: Particle position instances are set in ParticlePositionPool for re-use.
     * 
     * @param aSimulationStep Simulation step
     * @param aTemperature Temperature
     * @param anUpotDpd DPD potential energy
     * @param anUpotBond Bond potential energy
     * @param anUpotElectrostatics Electrostatics potential energy
     * @param anUpotTotal Total potential energy (= aUpotDdpd + aUpotBond + aUpotElectrostatics)
     * @param anUkin Kinetic energy
     * @param anUtotal Total energy
     * @param aSurfaceTensionAlongX Surface tension along x axis
     * @param aSurfaceTensionAlongY Surface tension along y axis
     * @param aSurfaceTensionAlongZ Surface tension along z axis
     * @param aSurfaceTensionNorm Norm (magnitude) of surface tension
     * @param aDpdSurfaceTensionAlongX DPD surface tension along x axis
     * @param aDpdSurfaceTensionAlongY DPD surface tension along y axis
     * @param aDpdSurfaceTensionAlongZ DPD surface tension along z axis
     * @param aDpdSurfaceTensionNorm Norm (magnitude) of DPD surface tension
     * @param aMoleculeRgValues Molecule Rg values (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap Base molecule-particle to nearest-neighbor molecule-particle frequency map (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap Base molecule-particle to nearest-neighbor particle frequency map (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap Base molecule-particle to nearest-neighbor molecule frequency map (may be null)
     * @param aBaseMoleculeToNearestNeighborMoleculeFrequencyMap Base molecule to nearest-neighbor molecule frequency map (may be null)
     * @param aBaseMoleculeToNearestNeighborMoleculeFrequencyMap Base molecule to nearest-neighbor molecule-tuple frequency map (may be null)
     * @param aParticlePositions Particle positions
     */
    @Override
    public void setSimulationStepInformation(
        int aSimulationStep, 
        double aTemperature, 
        double anUpotDpd,
        double anUpotBond,
        double anUpotElectrostatics,
        double anUpotTotal, 
        double anUkin, 
        double anUtotal, 
        double aSurfaceTensionAlongX,
        double aSurfaceTensionAlongY,
        double aSurfaceTensionAlongZ,
        double aSurfaceTensionNorm,
        double aDpdSurfaceTensionAlongX,
        double aDpdSurfaceTensionAlongY,
        double aDpdSurfaceTensionAlongZ,
        double aDpdSurfaceTensionNorm,
        MoleculeRgValue[] aMoleculeRgValues,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeToNearestNeighborMoleculeFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
        ParticlePosition[] aParticlePositions
    ) {
        String tmpSimulationStepString = String.valueOf(aSimulationStep);
        // Set temperature
        this.temperatureList.add(tmpSimulationStepString);
        this.temperatureList.add(String.valueOf(aTemperature));
        // Set DPD potential energy
        this.uPotDpdList.add(tmpSimulationStepString);
        this.uPotDpdList.add(String.valueOf(anUpotDpd));
        // Set bond potential energy
        this.uPotBondList.add(tmpSimulationStepString);
        this.uPotBondList.add(String.valueOf(anUpotBond));
        // Set electrostatics potential energy
        this.uPotElectrostaticsList.add(tmpSimulationStepString);
        this.uPotElectrostaticsList.add(String.valueOf(anUpotElectrostatics));
        // Set total potential energy
        this.uPotTotalList.add(tmpSimulationStepString);
        this.uPotTotalList.add(String.valueOf(anUpotTotal));
        // Set kinetic energy
        this.uKinList.add(tmpSimulationStepString);
        this.uKinList.add(String.valueOf(anUkin));
        // Set total energy
        this.uTotalList.add(tmpSimulationStepString);
        this.uTotalList.add(String.valueOf(anUtotal));
        // Set surface tension along x axis
        this.surfaceTensionAlongXList.add(tmpSimulationStepString);
        this.surfaceTensionAlongXList.add(String.valueOf(aSurfaceTensionAlongX));
        // Set surface tension along y axis
        this.surfaceTensionAlongYList.add(tmpSimulationStepString);
        this.surfaceTensionAlongYList.add(String.valueOf(aSurfaceTensionAlongY));
        // Set surface tension along z axis
        this.surfaceTensionAlongZList.add(tmpSimulationStepString);
        this.surfaceTensionAlongZList.add(String.valueOf(aSurfaceTensionAlongZ));
        // Set surface tension norm
        this.surfaceTensionNormList.add(tmpSimulationStepString);
        this.surfaceTensionNormList.add(String.valueOf(aSurfaceTensionNorm));
        // Set DPD surface tension along x axis
        this.dpdSurfaceTensionAlongXList.add(tmpSimulationStepString);
        this.dpdSurfaceTensionAlongXList.add(String.valueOf(aDpdSurfaceTensionAlongX));
        // Set DPD surface tension along y axis
        this.dpdSurfaceTensionAlongYList.add(tmpSimulationStepString);
        this.dpdSurfaceTensionAlongYList.add(String.valueOf(aDpdSurfaceTensionAlongY));
        // Set DPD surface tension along z axis
        this.dpdSurfaceTensionAlongZList.add(tmpSimulationStepString);
        this.dpdSurfaceTensionAlongZList.add(String.valueOf(aDpdSurfaceTensionAlongZ));
        // Set DPD surface tension norm
        this.dpdSurfaceTensionNormList.add(tmpSimulationStepString);
        this.dpdSurfaceTensionNormList.add(String.valueOf(aDpdSurfaceTensionNorm));
        // <editor-fold defaultstate="collapsed" desc="Evaluate Rg calculation results if necessary">
        if (aMoleculeRgValues != null) {
            // <editor-fold defaultstate="collapsed" desc="Initialize if necessary">
            if (this.rgMoleculeNames == null) {
                this.rgMoleculeNames = new String[aMoleculeRgValues.length];
                for (int i = 0; i < aMoleculeRgValues.length; i++) {
                    this.rgMoleculeNames[i] = aMoleculeRgValues[i].getMoleculeName();
                }
                this.rgValueLists = new LinkedList[this.rgMoleculeNames.length];
                for (int i = 0; i < this.rgValueLists.length; i++) {
                    this.rgValueLists[i] = this.getInitializedList(this.getRgFilePathname(this.rgMoleculeNames[i]));
                    // Important: Add molecule name first
                    this.rgValueLists[i].add(this.rgMoleculeNames[i]);
                }
            }
            // </editor-fold>
            for (int i = 0; i < this.rgMoleculeNames.length; i++) {
                this.rgValueLists[i].add(tmpSimulationStepString);
                this.rgValueLists[i].add(String.valueOf(aMoleculeRgValues[i].getRgValue()));
            }
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Set nearest-neighbor frequencies if necessary">
        if (aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap != null) {
            if (this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap == null) {
                this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap = 
                    this.getInitializedBaseToNearestNeighborMap(
                        this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname, 
                        aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap.size()
                    );
            }
            this.updateBaseToNearestNeighborMap(
                aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap,
                this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap,
                this.simulationStepList,
                tmpSimulationStepString
            );
        }
        if (aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap != null) {
            if (this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap == null) {
                this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap = 
                    this.getInitializedBaseToNearestNeighborMap(
                        this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname, 
                        aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap.size()
                    );
            }
            this.updateBaseToNearestNeighborMap(
                aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap,
                this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap,
                this.simulationStepList,
                tmpSimulationStepString
            );
        }
        if (aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap != null) {
            if (this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap == null) {
                this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap = 
                    this.getInitializedBaseToNearestNeighborMap(
                        this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname, 
                        aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap.size()
                    );
            }
            this.updateBaseToNearestNeighborMap(
                aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap,
                this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap,
                this.simulationStepList,
                tmpSimulationStepString
            );
        }
        if (aBaseMoleculeToNearestNeighborMoleculeFrequencyMap != null) {
            if (this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMap == null) {
                this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMap = 
                    this.getInitializedBaseToNearestNeighborMap(
                        this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname, 
                        aBaseMoleculeToNearestNeighborMoleculeFrequencyMap.size()
                    );
            }
            this.updateBaseToNearestNeighborMap(
                aBaseMoleculeToNearestNeighborMoleculeFrequencyMap,
                this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMap,
                this.simulationStepList,
                tmpSimulationStepString
            );
        }
        if (aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap != null) {
            if (this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap == null) {
                this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap = 
                    this.getInitializedBaseToNearestNeighborMap(
                        this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname, 
                        aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap.size()
                    );
            }
            this.updateBaseToNearestNeighborMap(
                aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
                this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap,
                this.simulationStepList,
                tmpSimulationStepString
            );
        }
        // </editor-fold>
        // Set simulation step AFTER nearest-neighbor settings
        this.simulationStepList.add(tmpSimulationStepString);
        // <editor-fold defaultstate="collapsed" desc="Set and write particle positions">
        if (this.simulationStepParticlePositionsDirectoryPath != null) {
            this.setCount++;
            String tmpFileEnding;
            tmpFileEnding = Strings.GZIP_FILE_ENDING;
            String tmpParticlePositionsSimulationStepFilePathname = 
                String.format(FileOutputStrings.PARTICLE_POSITIONS_SIMULATION_STEP_FILE_PREFIX_FORMAT, this.simulationStepParticlePositionsDirectoryPath + File.separatorChar, String.valueOf(aSimulationStep)) + tmpFileEnding;
            if (aSimulationStep > this.lastSimulationStep) {
                this.lastSimulationStep = aSimulationStep;
                this.lastSimulationStepParticlePositionsFilePathname = tmpParticlePositionsSimulationStepFilePathname;
            }
            this.executorService.execute(
                new ParticlePositionsWriteTask(
                    tmpParticlePositionsSimulationStepFilePathname, 
                    aParticlePositions, 
                    this.particlePositionsfilePathnameQueue, 
                    this.particlePositionPool,
                    this.numberOfAfterDecimalDigitsForParticlePositions
                )
            );
        }
        // </editor-fold>
    }
    
    /**
     * Sets restart info
     * 
     * @param aRestartInfo Restart info
     */
    @Override
    public void setRestartInfo(RestartInfo aRestartInfo) {
        this.restartInfo = aRestartInfo;
    }

    /**
     * Particle position pool
     * 
     * @param aParticlePositionPool Particle position pool
     */
    @Override
    public void setParticlePositionPool(ParticlePositionPool aParticlePositionPool) {
        this.particlePositionPool = aParticlePositionPool;
    }
    
    /**
     * Finishes output
     */
    @Override
    public void finish() {
        this.writeSimulationStepProperties();
        if (this.executorService != null) {
            this.executorService.shutdown();
        }
        if (this.setCount != this.particlePositionsfilePathnameQueue.size()) {
            try {
                if (this.executorService != null) {
                    // 1 minute for time-out is expected to be sufficient for all practical purposes
                    this.executorService.awaitTermination(1L, TimeUnit.MINUTES);
                }
            } catch (Exception anException) {
                return;
            }
        }
        // Copy last step file to final file if necessary
        if (this.lastSimulationStepParticlePositionsFilePathname != null) {
            String tmpFileEnding;
            tmpFileEnding = Strings.GZIP_FILE_ENDING;
            Utils.copySingleFile(
                this.lastSimulationStepParticlePositionsFilePathname, 
                this.outputDirectoryPath + File.separatorChar + FileOutputStrings.PARTICLE_POSITIONS_FINAL_FILENAME_PREFIX + tmpFileEnding
            );
        }
        // Write restart info if necessary
        if (this.restartInfo != null) {
            String tmpFileEnding;
            tmpFileEnding = Strings.GZIP_FILE_ENDING;
            this.restartInfo.writeToFile(this.outputDirectoryPath + File.separatorChar + FileOutputStrings.RESTART_INFO_FILENAME_PREFIX + tmpFileEnding);
        }
    }

    /**
     * Writes simulation step properties
     */
    public void writeSimulationStepProperties() {
        Utils.deleteSingleFile(this.simulationStepFilePathname);
        Utils.writeStringArrayToFile(this.simulationStepList.toArray(new String[0]), this.simulationStepFilePathname);

        Utils.deleteSingleFile(this.temperatureFilePathname);
        Utils.writeStringArrayToFile(this.temperatureList.toArray(new String[0]), this.temperatureFilePathname);
        
        Utils.deleteSingleFile(this.uPotDpdFilePathname);
        Utils.writeStringArrayToFile(this.uPotDpdList.toArray(new String[0]), this.uPotDpdFilePathname);
        Utils.deleteSingleFile(this.uPotBondFilePathname);
        Utils.writeStringArrayToFile(this.uPotBondList.toArray(new String[0]), this.uPotBondFilePathname);
        Utils.deleteSingleFile(this.uPotElectrostaticsFilePathname);
        Utils.writeStringArrayToFile(this.uPotElectrostaticsList.toArray(new String[0]), this.uPotElectrostaticsFilePathname);
        Utils.deleteSingleFile(this.uPotTotalFilePathname);
        Utils.writeStringArrayToFile(this.uPotTotalList.toArray(new String[0]), this.uPotTotalFilePathname);
        Utils.deleteSingleFile(this.uKinFilePathname);
        Utils.writeStringArrayToFile(this.uKinList.toArray(new String[0]), this.uKinFilePathname);
        Utils.deleteSingleFile(this.uTotalFilePathname);
        Utils.writeStringArrayToFile(this.uTotalList.toArray(new String[0]), this.uTotalFilePathname);

        Utils.deleteSingleFile(this.surfaceTensionAlongXFilePathname);
        Utils.writeStringArrayToFile(this.surfaceTensionAlongXList.toArray(new String[0]), this.surfaceTensionAlongXFilePathname);
        Utils.deleteSingleFile(this.surfaceTensionAlongYFilePathname);
        Utils.writeStringArrayToFile(this.surfaceTensionAlongYList.toArray(new String[0]), this.surfaceTensionAlongYFilePathname);
        Utils.deleteSingleFile(this.surfaceTensionAlongZFilePathname);
        Utils.writeStringArrayToFile(this.surfaceTensionAlongZList.toArray(new String[0]), this.surfaceTensionAlongZFilePathname);
        Utils.deleteSingleFile(this.surfaceTensionNormFilePathname);
        Utils.writeStringArrayToFile(this.surfaceTensionNormList.toArray(new String[0]), this.surfaceTensionNormFilePathname);

        Utils.deleteSingleFile(this.dpdSurfaceTensionAlongXFilePathname);
        Utils.writeStringArrayToFile(this.dpdSurfaceTensionAlongXList.toArray(new String[0]), this.dpdSurfaceTensionAlongXFilePathname);
        Utils.deleteSingleFile(this.dpdSurfaceTensionAlongYFilePathname);
        Utils.writeStringArrayToFile(this.dpdSurfaceTensionAlongYList.toArray(new String[0]), this.dpdSurfaceTensionAlongYFilePathname);
        Utils.deleteSingleFile(this.dpdSurfaceTensionAlongZFilePathname);
        Utils.writeStringArrayToFile(this.dpdSurfaceTensionAlongZList.toArray(new String[0]), this.dpdSurfaceTensionAlongZFilePathname);
        Utils.deleteSingleFile(this.dpdSurfaceTensionNormFilePathname);
        Utils.writeStringArrayToFile(this.dpdSurfaceTensionNormList.toArray(new String[0]), this.dpdSurfaceTensionNormFilePathname);
        
        if (this.rgMoleculeNames != null) {
            for (int i = 0; i < this.rgMoleculeNames.length; i++) {
                String tmpRgFilePathname = this.getRgFilePathname(this.rgMoleculeNames[i]);
                Utils.deleteSingleFile(tmpRgFilePathname);
                Utils.writeStringArrayToFile(this.rgValueLists[i].toArray(new String[0]), tmpRgFilePathname);
            }
        }
        
        Utils.deleteSingleFile(this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname);
        Utils.writeBaseToNearestNeighborStepFrequencyMapToFile(
            this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMap, 
            this.baseMoleculeParticleToNearestNeighborMoleculeParticleStepFrequencyMapFilePathname
        );
        Utils.deleteSingleFile(this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname);
        Utils.writeBaseToNearestNeighborStepFrequencyMapToFile(
            this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMap, 
            this.baseMoleculeParticleToNearestNeighborParticleStepFrequencyMapFilePathname
        );
        Utils.deleteSingleFile(this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname);
        Utils.writeBaseToNearestNeighborStepFrequencyMapToFile(
            this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMap, 
            this.baseMoleculeParticleToNearestNeighborMoleculeStepFrequencyMapFilePathname
        );
        Utils.deleteSingleFile(this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname);
        Utils.writeBaseToNearestNeighborStepFrequencyMapToFile(
            this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMap, 
            this.baseMoleculeToNearestNeighborMoleculeStepFrequencyMapFilePathname
        );
        Utils.deleteSingleFile(this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname);
        Utils.writeBaseToNearestNeighborStepFrequencyMapToFile(
            this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMap, 
            this.baseMoleculeToNearestNeighborMoleculeTupleStepFrequencyMapFilePathname
        );
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Queue with particle positions file pathnames
     * 
     * @return Queue with particle positions file pathnames (may be empty)
     */
    public ConcurrentLinkedQueue<String> getParticlePositionsFilePathnameQueue() {
        return this.particlePositionsfilePathnameQueue;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Returns initialised list (see code)
     * 
     * @param aFilePathname File pathname of previous list data (may be null or unavailable)
     * @return Initialised list
     */
    private LinkedList<String> getInitializedList(String aFilePathname) {
        if (aFilePathname != null && !aFilePathname.isEmpty() && (new File(aFilePathname)).isFile()) {
            String[] tmpStringArray = Utils.readStringArrayFromFile(aFilePathname);
            if (tmpStringArray[0].equals(Strings.VERSION_1_0_0)) {
                LinkedList<String> tmpStringList = new LinkedList<>();
                for (int i = 0; i < tmpStringArray.length; i++) {
                    tmpStringList.add(tmpStringArray[i]);
                }
                return tmpStringList;
            }
        }
        LinkedList<String> tmpStringList = new LinkedList<>();
        tmpStringList.add(Strings.VERSION_1_0_0);
        return tmpStringList;
    }
    
    /**
     * Pathname of Rg file
     * NOTE: NO checks are performed.
     * 
     * @param aMoleculeName Molecule name
     * @return Pathname of Rg file
     */
    private String getRgFilePathname(String aMoleculeName) {
        return this.radiusOfGyrationDirectoryPath + File.separatorChar + FileOutputStrings.RADIUS_OF_GYRATION_FILENAME_PREFIX + aMoleculeName + FileOutputStrings.TEXT_FILE_ENDING;
    }
    
    /**
     * Updates base to nearest-neighbor map (see code)
     * NOTE: NO checks are performed.
     * 
     * @param aBaseToNearestNeighborFrequencyMap Base to nearest-neighbor frequency map (NOT allowed to be null, is NOT changed)
     * @param aBaseToNearestNeighborStepFrequencyMap Base to nearest-neighbor step-frequency map (NOT allowed to be null, may be changed)
     * @param aSimulationStepList Simulation step list WIHTOUT aCurrentSimulationStepString (NOT allowed to be null, is NOT changed)
     * @param aCurrentSimulationStepString Current simulation step string (NOT allowed to be null/empty)
     */
    private void updateBaseToNearestNeighborMap(
        HashMap<String, HashMap<String, Integer>> aBaseToNearestNeighborFrequencyMap,
        HashMap<String, HashMap<String, LinkedList<String>>> aBaseToNearestNeighborStepFrequencyMap,
        LinkedList<String> aSimulationStepList,
        String aCurrentSimulationStepString
    )
    {
        String tmpZeroString = String.valueOf(0);
        for (String tmpBase : aBaseToNearestNeighborFrequencyMap.keySet()) {
            HashMap<String, Integer> tmpNearestNeighborToFrequencyMap = aBaseToNearestNeighborFrequencyMap.get(tmpBase);
            if (aBaseToNearestNeighborStepFrequencyMap.containsKey(tmpBase)) {
                // <editor-fold defaultstate="collapsed" desc="aBaseToNearestNeighborStepFrequencyMap contains tmpBase">
                HashMap<String, LinkedList<String>> tmpNearestNeighborToStepFrequencyMap = aBaseToNearestNeighborStepFrequencyMap.get(tmpBase);
                for (String tmpNearestNeighbor : tmpNearestNeighborToFrequencyMap.keySet()) {
                    if (tmpNearestNeighborToStepFrequencyMap.containsKey(tmpNearestNeighbor)) {
                        LinkedList<String> tmpStepFrequencyList = tmpNearestNeighborToStepFrequencyMap.get(tmpNearestNeighbor);
                        tmpStepFrequencyList.add(aCurrentSimulationStepString);
                        tmpStepFrequencyList.add(String.valueOf(tmpNearestNeighborToFrequencyMap.get(tmpNearestNeighbor)));
                    } else {
                        LinkedList<String> tmpStepFrequencyList = this.getStepFrequencyList(aSimulationStepList);
                        tmpStepFrequencyList.add(aCurrentSimulationStepString);
                        tmpStepFrequencyList.add(String.valueOf(tmpNearestNeighborToFrequencyMap.get(tmpNearestNeighbor)));
                        tmpNearestNeighborToStepFrequencyMap.put(tmpNearestNeighbor, tmpStepFrequencyList);
                    }
                }
                for (String tmpExistingNearestNeighbor : tmpNearestNeighborToStepFrequencyMap.keySet()) {
                    if (!tmpNearestNeighborToFrequencyMap.containsKey(tmpExistingNearestNeighbor)) {
                        LinkedList<String> tmpStepFrequencyList = tmpNearestNeighborToStepFrequencyMap.get(tmpExistingNearestNeighbor);
                        tmpStepFrequencyList.add(aCurrentSimulationStepString);
                        tmpStepFrequencyList.add(tmpZeroString);
                    }
                }
                // </editor-fold>
            } else {
                // <editor-fold defaultstate="collapsed" desc="aBaseToNearestNeighborStepFrequencyMap does NOT contain tmpBase">
                HashMap<String, LinkedList<String>> tmpNearestNeighborToStepFrequencyMap = new HashMap<>(tmpNearestNeighborToFrequencyMap.size());
                for (String tmpNearestNeighbor : tmpNearestNeighborToFrequencyMap.keySet()) {
                    LinkedList<String> tmpStepFrequencyList = this.getStepFrequencyList(aSimulationStepList);
                    tmpStepFrequencyList.add(aCurrentSimulationStepString);
                    tmpStepFrequencyList.add(String.valueOf(tmpNearestNeighborToFrequencyMap.get(tmpNearestNeighbor)));
                    tmpNearestNeighborToStepFrequencyMap.put(tmpNearestNeighbor, tmpStepFrequencyList);
                }
                aBaseToNearestNeighborStepFrequencyMap.put(tmpBase, tmpNearestNeighborToStepFrequencyMap);
                // </editor-fold>
            }
        }
    }

    /**
     * Support method for method updateBaseToNearestNeighborMap().
     * NOTE: NO checks are performed.
     * 
     * @param aSimulationStepList Simulation step list (NOT allowed to be null, is NOT changed)
     */
    private LinkedList<String> getStepFrequencyList(LinkedList<String> aSimulationStepList) {
        LinkedList<String> tmpStepFrequencyList = new LinkedList<>();
        if (aSimulationStepList.size() > 0) {
            String tmpZeroString = String.valueOf(0);
            // NOTE: FIRST line is version
            boolean tmpIsFirstLine = true;
            for (String tmpElement : aSimulationStepList) {
                if (tmpIsFirstLine) {
                    tmpStepFrequencyList.add(tmpElement);
                    tmpIsFirstLine = false;
                } else {
                    tmpStepFrequencyList.add(tmpElement);
                    tmpStepFrequencyList.add(tmpZeroString);
                }
            }
        }
        return tmpStepFrequencyList;
    }

    /**
     * Returns initialised BaseToNearestNeighborMap, see code.
     * 
     * @param aFilePathname File pathname of BaseToNearestNeighborMap
     * @param aSize Size of BaseToNearestNeighborMap hash map 
     */
    private HashMap<String, HashMap<String, LinkedList<String>>> getInitializedBaseToNearestNeighborMap(String aFilePathname, int aSize) {
        if (aFilePathname != null && !aFilePathname.isEmpty() && (new File(aFilePathname)).isFile()) {
            return Utils.readBaseToNearestNeighborStepFrequencyMapFromFile(aFilePathname);
        } else {
            return new HashMap<>(aSize);
        }
    }
    // </editor-fold>
    
}
