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
package de.gnwi.jdpd.samples.interactions;

import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.utilities.CellBox;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.interfaces.IParticlePairInteractionCalculator;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.utilities.CalculatorUtils;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.RandomAdderGroup;
import de.gnwi.jdpd.utilities.Utils;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * Particle pair interaction calculator that uses CellBox and the cell linked-list 
 * method to determine particle pair interactions that are restricted to a defined 
 * cut-off length.
 * 
 * For basic concept and discussion see
 * M.P. Allen and D.J. Tildesley
 * Computer Simulation of Liquids
 * Clarendon Press
 * Oxford 1987
 * Chapter 5.3
 * 
 * The cell chunk parallelisation implementation allows the use of fast NON 
 * thread-safe objects and accumulation items like non thread-safe random number 
 * generators, non thread-safe double adders and standard (NON atomic) double 
 * arrays.
 * 
 * @author Achim Zielesny
 */
public abstract class ParticlePairInteractionCalculator extends CellBox implements IParticlePairInteractionCalculator {

    // <editor-fold defaultstate="collapsed" desc="Public Enums">
    /**
     * Calculation mode
     */
    public enum CellBasedCalculationMode {
            
        /**
         * Particle cell assignments are calculated
         */
        WITH_PARTICLE_CELL_ASSIGNMENTS,
        /**
         * Particle cell assignments are NOT calculated, current assignments are used
         */
        WITHOUT_PARTICLE_CELL_ASSIGNMENTS,
        /**
         * Cache is used
         */
        WITH_CACHE,
        /**
         * Undefined
         */
        UNDEFINED

    }  
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class CellChunkCalculationTask">
    /**
     * Callable for calculation of particle pair interactions of cell chunks in a 
     * parallelised manner with executor service
     */
    private class CellChunkCalculationTask implements Callable<Boolean> {
        
        // <editor-fold defaultstate="collapsed" desc="Private final class variables">
        /**
         * Simulation logger
         */
        private final ILogger simulationLogger;
        
        /**
         * Cell chunks that are parallelisation-safe
         */
        private final int[] parallelizationSafeCellChunk;

        /**
         * Start index for parallelizationSafeCellChunk
         */
        private final int startIndex;

        /**
         * End index for parallelizationSafeCellChunk
         */
        private final int endIndex;
        
        /**
         * Particle x-positions
         */
        private final double[] r_x;
        
        /**
         * Particle y-positions
         */
        private final double[] r_y;
        
        /**
         * Particle z-positions
         */
        private final double[] r_z;
        
        /**
         * Parameters
         */
        private final Parameters parameters;
        
        /**
         * Random adder group
         */
        private final RandomAdderGroup randomAdderGroup;
        
        /**
         * ParticlePairDistanceParameters instance for cache
         */
        private final ParticlePairDistanceParameters particlePairDistanceParameters;
        
        /**
         * Particle pair interaction calculator
         */
        private final ParticlePairInteractionCalculator particlePairInteractionCalculator;
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Constructor">
        /**
         * Constructor
         * 
         * @param aSimulationLogger Simulation logger
         * @param aParallelizationSafeCellChunk Cell chunks that are parallelisation-safe
         * @param aStartIndex Start index for parallelizationSafeCellChunk
         * @param anEndIndex End index for parallelizationSafeCellChunk
         * @param aR_x Current x-components of particle positions in simulation box
         * @param aR_y Current y-components of particle positions in simulation box
         * @param aR_z Current z-components of particle positions in simulation box
         * @param aRandomAdderGroup Random adder group
         * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance for cache
         * @param aParameters Parameters (may be null)
         * @param aParticlePairInteractionCalculator Particle pair interaction calculator
         */
        public CellChunkCalculationTask(
            ILogger aSimulationLogger,
            int[] aParallelizationSafeCellChunk,
            int aStartIndex, 
            int anEndIndex,
            double[] aR_x,
            double[] aR_y,
            double[] aR_z,
            RandomAdderGroup aRandomAdderGroup,
            ParticlePairDistanceParameters aParticlePairDistanceParameters,
            Parameters aParameters,
            ParticlePairInteractionCalculator aParticlePairInteractionCalculator
        ) {
            this.simulationLogger = aSimulationLogger;
            this.parallelizationSafeCellChunk = aParallelizationSafeCellChunk;
            this.startIndex = aStartIndex;
            this.endIndex = anEndIndex;
            this.r_x = aR_x;
            this.r_y = aR_y;
            this.r_z = aR_z;
            this.randomAdderGroup = aRandomAdderGroup;
            this.particlePairDistanceParameters = aParticlePairDistanceParameters;
            this.parameters = aParameters;
            this.particlePairInteractionCalculator = aParticlePairInteractionCalculator;
        }
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Public overridden methods">
        /**
         * Overridden call method
         * 
         * @return Always true
         */
        @Override
        public Boolean call() {
            try {
                for (int i = this.startIndex; i < this.endIndex; i++) {
                    int tmpCellIndex = this.parallelizationSafeCellChunk[i];
                    this.particlePairInteractionCalculator.calculateSingleCellParticlePairInteractions(
                        tmpCellIndex,
                        this.r_x,
                        this.r_y,
                        this.r_z,
                        this.randomAdderGroup,
                        this.particlePairDistanceParameters,
                        this.parameters
                    );
                }
                return true;
            } catch (Exception anException) {
                // <editor-fold defaultstate="collapsed" desc="Exception logging">
                this.simulationLogger.appendException("CellChunkCalculationTask.call", Utils.getStacktrace(anException));
                // </editor-fold>
                throw anException;
            }
        }
        // </editor-fold>
        
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class CacheCalculationTask">
    /**
     * Callable for calculation of particle pair interactions in cache in a 
     * parallelised manner with executor service
     */
    private class CacheCalculationTask implements Callable<Boolean> {
        
        // <editor-fold defaultstate="collapsed" desc="Private final class variables">
        /**
         * Simulation logger
         */
        private final ILogger simulationLogger;

        /**
         * Parameters
         */
        private final Parameters parameters;
        
        /**
         * Random adder group
         */
        private final RandomAdderGroup randomAdderGroup;
        
        /**
         * ParticlePairDistanceParameters instance for cache
         */
        private final ParticlePairDistanceParameters particlePairDistanceParameters;
        
        /**
         * Particle pair interaction calculator
         */
        private final ParticlePairInteractionCalculator particlePairInteractionCalculator;
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Constructor">
        /**
         * Constructor
         * 
         * @param aSimulationLogger Simulation logger
         * @param aParallelizationSafeCellChunk Cell chunks that are parallelisation-safe
         * @param aStartIndex Start index for parallelizationSafeCellChunk
         * @param anEndIndex End index for parallelizationSafeCellChunk
         * @param aR_x Current x-components of particle positions in simulation box
         * @param aR_y Current y-components of particle positions in simulation box
         * @param aR_z Current z-components of particle positions in simulation box
         * @param aRandomAdderGroup Random adder group
         * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance for cache
         * @param aParameters Parameters (may be null)
         * @param aParticlePairInteractionCalculator Particle pair interaction calculator
         */
        public CacheCalculationTask(
            ILogger aSimulationLogger,
            RandomAdderGroup aRandomAdderGroup,
            ParticlePairDistanceParameters aParticlePairDistanceParameters,
            Parameters aParameters,
            ParticlePairInteractionCalculator aParticlePairInteractionCalculator
        ) {
            this.simulationLogger = aSimulationLogger;
            this.randomAdderGroup = aRandomAdderGroup;
            this.particlePairDistanceParameters = aParticlePairDistanceParameters;
            this.parameters = aParameters;
            this.particlePairInteractionCalculator = aParticlePairInteractionCalculator;
        }
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Public overridden methods">
        /**
         * Overridden call method
         * 
         * @return Always true
         */
        @Override
        public Boolean call() {
            try{
                int[] tmpParticleIndices1 = this.particlePairDistanceParameters.getParticleIndices_i();
                int[] tmpParticleIndices2 = this.particlePairDistanceParameters.getParticleIndices_j();
                double[] tmpRij_x = this.particlePairDistanceParameters.getArray_Rij_x();
                double[] tmpRij_y = this.particlePairDistanceParameters.getArray_Rij_y();
                double[] tmpRij_z = this.particlePairDistanceParameters.getArray_Rij_z();
                double[] tmpRij_Square = this.particlePairDistanceParameters.getArray_Rij_Square();
                double[] tmpRij = this.particlePairDistanceParameters.getArray_Rij();
                for (int i = 0; i < this.particlePairDistanceParameters.getSize(); i++) {
                    this.particlePairInteractionCalculator.calculateParticlePairInteraction(
                        tmpParticleIndices1[i],
                        tmpParticleIndices2[i], 
                        tmpRij_x[i],
                        tmpRij_y[i],
                        tmpRij_z[i],
                        tmpRij_Square[i],
                        tmpRij[i],
                        this.randomAdderGroup,
                        this.parameters
                    );
                }
                return true;
            } catch (Exception anException) {
                // <editor-fold defaultstate="collapsed" desc="Exception logging">
                this.simulationLogger.appendException("CacheCalculationTask.call", Utils.getStacktrace(anException));
                // </editor-fold>
                throw anException;
            }
        }
        // </editor-fold>
        
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Factory for new objects
     */
    private final Factory factory;
    
    /**
     * Square of cut-off length
     */
    private final double cutOffLengthSquare;
    
    /**
     * Parallelisation info
     */
    private final ParallelizationInfo parallelizationInfo;

    /**
     * Random number seed
     */
    private final AtomicInteger randomNumberSeed;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Head particle index of cell
     */
    private int[] headParticleIndexOfCell;
    
    /**
     * Next particle index of cell
     */
    private int[] nextParticleIndexOfCell;
    
    /**
     * Cache for ParticlePairDistanceParameterss
     */
    private ParticlePairDistanceParameters[][] particlePairDistanceParametersCache;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected final class variables">
    /**
     * Calculator utils
     */
    protected final CalculatorUtils calculatorUtils;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected class variables">
    /**
     * True: Cache for ParticlePairDistanceParameterss is active, false: Otherwise
     */
    protected boolean isParticlePairDistanceParametersCacheActive;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aFactory Factory for new objects
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected ParticlePairInteractionCalculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException 
    {
        this(
            aFactory,
            aSimulationLogger,
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength,
            null,
            null,
            aParallelizationInfo,
            aRandomNumberSeed
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairInteractionCalculator.Constructor WITHOUT aCellNeigbours, aParallelizationSafeCellChunks");
        // </editor-fold>
    }

    /**
     * Constructor
     * 
     * @param aFactory Factory for new objects
     * @param aSimulationLogger Simulation logger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aCellNeigbours Cell neighbours (may be null then cell neighbours are determined)
     * @param aParallelizationSafeCellChunks Cell chunks that are parallelisation-safe (may be null then cell chunks are determined)
     * @param aParallelizationInfo Parallelisation info
     * @param aRandomNumberSeed Random number seed
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected ParticlePairInteractionCalculator(
        Factory aFactory,
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        int[][] aCellNeigbours,
        int[][] aParallelizationSafeCellChunks,
        ParallelizationInfo aParallelizationInfo,
        AtomicInteger aRandomNumberSeed
    ) throws IllegalArgumentException 
    {
        super(
            aSimulationLogger,
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength,
            aCellNeigbours,
            aParallelizationSafeCellChunks
        );
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aFactory == null) {
            throw new IllegalArgumentException("ParticlePairInteractionCalculator.Constructor: aFactory is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("ParticlePairInteractionCalculator.Constructor: aParallelizationInfo is null.");
        }
        // </editor-fold>
        try {
            this.factory = aFactory;
            this.parallelizationInfo = aParallelizationInfo;
            this.randomNumberSeed = aRandomNumberSeed;
            this.calculatorUtils = new CalculatorUtils();
            this.cutOffLengthSquare = this.cutOffLength * this.cutOffLength;
            this.nextParticleIndexOfCell = null;
            this.headParticleIndexOfCell = new int[this.cellNumber];
            this.isParticlePairDistanceParametersCacheActive = false;
            this.particlePairDistanceParametersCache = null;
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            this.simulationLogger.appendMethodCall("ParticlePairInteractionCalculator.Constructor: FULL");
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.Constructor", Utils.getStacktrace(anException));
            // </editor-fold>
            throw anException;
        }
    }
    
    /**
     * Constructor that clones aParticlePairInteractionCalculator (i.e. uses the same
 getBoxSize(), getPeriodicBoundaries() and getCellNeigbours() etc.)
 NOTE: getHeadParticleIndexOfCellArray() and getNextParticleIndexOfCellArray are
       NOT cloned but may be set manually with method 
       this.setParticleCellAssignments() if available.
     * 
     * @param aParticlePairInteractionCalculator ParticlePairInteractionCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected ParticlePairInteractionCalculator(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) throws IllegalArgumentException {
        this(
            aParticlePairInteractionCalculator.getFactory(),
            aParticlePairInteractionCalculator.getSimulationLogger(), 
            aParticlePairInteractionCalculator.getBoxSize(), 
            aParticlePairInteractionCalculator.getPeriodicBoundaries(), 
            aParticlePairInteractionCalculator.getCutOffLength(),
            aParticlePairInteractionCalculator.getCellNeigbours(),
            aParticlePairInteractionCalculator.getParallelizationSafeCellChunks(),
            aParticlePairInteractionCalculator.getParallelizationInfo(),
            aParticlePairInteractionCalculator.getRandomNumberSeed()
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("ParticlePairInteractionCalculator.Constructor WITH aParticlePairInteractionCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    // <editor-fold defaultstate="collapsed" desc="- Calculation method">
    /**
     * Calculates particle pair interactions with a distance smaller than 
     * cut-off length based on calculation mode.
     * (No checks are performed)
     * NOTE: This method scales linearly with the number of particles.
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters
     * @param aCalculationMode Calculation mode
     * @return True: Operation successful, false: Otherwise
     * @throws IllegalArgumentException Thrown if calculation mode is illegal
     */
    @Override
    public boolean calculateParticlePairInteractions(
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters,
        CellBasedCalculationMode aCalculationMode
    ) {
        switch(aCalculationMode) {
            case WITH_PARTICLE_CELL_ASSIGNMENTS:
                return calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
                    aR_x,
                    aR_y,
                    aR_z,
                    aParameters
                );
            case WITHOUT_PARTICLE_CELL_ASSIGNMENTS:
                return calculateCellBasedParticlePairInteractionsWithoutParticleCellAssignments(
                    aR_x,
                    aR_y,
                    aR_z,
                    aParameters
                );
            case WITH_CACHE:
                return calculateCacheBasedParticlePairInteractions(
                    aParameters
                );
            case UNDEFINED:
                throw new IllegalArgumentException("ParticlePairInteractionCalculator.calculateCellBasedParticlePairInteractions: Calculation mode is undefined.");
            default:
                throw new IllegalArgumentException("ParticlePairInteractionCalculator.calculateCellBasedParticlePairInteractions: Unknown calculation mode.");
        }
    }    
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Cell-based calculation methods">    
    /**
     * Calculates particle pair interactions with a distance smaller than 
     * cut-off length based on cell partitioning. Particle cell assignments are 
     * calculated.
     * (No checks are performed)
     * NOTE: This method scales linearly with the number of particles.
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    public boolean calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters
    ) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Choose calculation method">
            if (!this.isCellBasedCalculation) {
                return this.calculateLoopBasedParticlePairInteractions(
                    aR_x,
                    aR_y,
                    aR_z, 
                    aParameters);
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Reset double adders">
            this.calculatorUtils.resetAdders();
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Initialise cell linked-list related particle index arrays">
            int tmpParticleNumber = aR_x.length;
            if (this.nextParticleIndexOfCell == null) {
                this.nextParticleIndexOfCell = new int[tmpParticleNumber];
            } else if (this.nextParticleIndexOfCell.length != tmpParticleNumber) {
                this.nextParticleIndexOfCell = new int[tmpParticleNumber];
            }
            Arrays.fill(this.nextParticleIndexOfCell, -1);
            // NOTE: this.headParticleIndexOfCell is already allocated by this.initialise()
            Arrays.fill(this.headParticleIndexOfCell, -1);
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Assign particle indices to cells">
            for (int tmpParticleIndex = 0; tmpParticleIndex < tmpParticleNumber; tmpParticleIndex++) {
                int tmpIndexX = this.getIndexAlongX(aR_x[tmpParticleIndex]);
                // NOTE: Correct for position exactly on right boundary
                if (tmpIndexX == this.cellNumberX) tmpIndexX = this.cellNumberX - 1;

                int tmpIndexY = this.getIndexAlongY(aR_y[tmpParticleIndex]);
                // NOTE: Correct for position exactly on right boundary
                if (tmpIndexY == this.cellNumberY) tmpIndexY = this.cellNumberY - 1;

                int tmpIndexZ = this.getIndexAlongZ(aR_z[tmpParticleIndex]);
                // NOTE: Correct for position exactly on right boundary
                if (tmpIndexZ == this.cellNumberZ) tmpIndexZ = this.cellNumberZ - 1;

                int tmpCellIndex = this.getCellIndex(tmpIndexX, tmpIndexY, tmpIndexZ);
                this.nextParticleIndexOfCell[tmpParticleIndex] = this.headParticleIndexOfCell[tmpCellIndex];
                this.headParticleIndexOfCell[tmpCellIndex] = tmpParticleIndex;
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Determine particle index pairs">
            return this.determineParticleIndexPairs(
                aR_x,
                aR_y,
                aR_z, 
                aParameters
            );
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
    }
    
    /**
     * Calculates particle pair interactions with a distance smaller than 
     * cut-off length based on cell partitioning.
     * NOTE: Particles MUST already be assigned to cells.
     * (No checks are performed)
     * NOTE: This method scales linearly with the number of particles.
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    public boolean calculateCellBasedParticlePairInteractionsWithoutParticleCellAssignments(
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters
    ) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Choose calculation method">
            if (!this.isCellBasedCalculation) {
                return this.calculateLoopBasedParticlePairInteractions(
                    aR_x,
                    aR_y,
                    aR_z, 
                    aParameters);
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Reset double adders">
            this.calculatorUtils.resetAdders();
            // </editor-fold>
            // NO initialisation of cell linked-list related particle index arrays
            // NO assignment of particle indices to cells
            // <editor-fold defaultstate="collapsed" desc="Determine particle index pairs">
            return this.determineParticleIndexPairs(
                aR_x,
                aR_y,
                aR_z, 
                aParameters
            );
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.calculateCellBasedParticlePairInteractionsWithoutParticleCellAssignments", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Cache-based calculation method">
    /**
     * Calculates particle pair interactions with a distance smaller than 
     * cut-off length based on cache.
     * NOTE: Cache MUST already be set.
     * (No checks are performed)
     * 
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    public boolean calculateCacheBasedParticlePairInteractions(Parameters aParameters) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Reset double adders">
            this.calculatorUtils.resetAdders();
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Calculate cache based particle index pairs">
            if (this.calculatorUtils.getAdderGroups() == null) {
                RandomAdderGroup[] tmpRandomAdderGroups = new RandomAdderGroup[this.particlePairDistanceParametersCache[0].length];
                for (int i = 0; i < this.particlePairDistanceParametersCache[0].length; i++) {
                    tmpRandomAdderGroups[i] = new RandomAdderGroup(this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet()));
                }
                this.calculatorUtils.setAdderGroups(tmpRandomAdderGroups);
            }
            if (this.calculatorUtils.getExecutorService() == null) {
                this.calculatorUtils.setExecutorService(Executors.newFixedThreadPool(this.particlePairDistanceParametersCache[0].length));
            }
            RandomAdderGroup[] tmpRandomAdderGroups = ((RandomAdderGroup[]) this.calculatorUtils.getAdderGroups());
            for (ParticlePairDistanceParameters[] tmpParticlePairDistanceParametersCachePerCellChunk : this.particlePairDistanceParametersCache) {
                LinkedList<CacheCalculationTask> tmpTaskList = new LinkedList<>();
                for (int k = 0; k < tmpParticlePairDistanceParametersCachePerCellChunk.length; k++) {
                    CacheCalculationTask newTask = 
                        new CacheCalculationTask(
                            this.simulationLogger,
                            tmpRandomAdderGroups[k], 
                            tmpParticlePairDistanceParametersCachePerCellChunk[k], 
                            aParameters, 
                            this
                        );
                    tmpTaskList.add(newTask);
                }
                this.calculatorUtils.getExecutorService().invokeAll(tmpTaskList);
            }
            return true;
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.calculateCacheBasedParticlePairInteractions", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Loop-based calculation methods">
    /**
     * Calculates particle pair interactions with a distance smaller than cut-off 
     * length by loop over all particle index pairs.
     * NOTE: No checks are performed
     * NOTE: This method scales quadratically with the number of particles.
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    public boolean calculateLoopBasedParticlePairInteractions(
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Loop over all particle index pairs">
            if (this.calculatorUtils.getAdderGroups() == null) {
                this.calculatorUtils.setAdderGroups(new RandomAdderGroup[] {new RandomAdderGroup(this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet()))});
            }
            this.calculatorUtils.resetAdders();
            RandomAdderGroup tmpRandomAdderGroup = ((RandomAdderGroup[]) this.calculatorUtils.getAdderGroups())[0];
            ParticlePairDistanceParameters tmpParticlePairDistanceParameters = null;
            if (this.isParticlePairDistanceParametersCacheActive && this.particlePairDistanceParametersCache == null) {
                this.particlePairDistanceParametersCache = new ParticlePairDistanceParameters[1][1];
                int tmpCapacity = aR_x.length;
                this.particlePairDistanceParametersCache[0][0] = new ParticlePairDistanceParameters(tmpCapacity);
            }
            if (this.isParticlePairDistanceParametersCacheActive) {
                tmpParticlePairDistanceParameters = this.particlePairDistanceParametersCache[0][0];
            }
            int tmpParticleNumber = aR_x.length;
            for (int tmpParticleIndex1 = 1; tmpParticleIndex1 < tmpParticleNumber; tmpParticleIndex1++) {
                for (int tmpParticleIndex2 = 0; tmpParticleIndex2 < tmpParticleIndex1; tmpParticleIndex2++) {
                    this.evaluateParticlePairInteraction(
                        tmpParticleIndex1, 
                        tmpParticleIndex2, 
                        aR_x, 
                        aR_y, 
                        aR_z, 
                        tmpRandomAdderGroup,
                        tmpParticlePairDistanceParameters,
                        aParameters
                    );
                }
            }
            // </editor-fold>
            return true;
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.calculateLoopBasedParticlePairInteractions", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Cache related methods">
    /**
     * Sets activity of cache for ParticlePairDistanceParameterss
     * 
     * @param aValue True: Cache for ParticlePairDistanceParameterss is active, false: Otherwise
     */
    @Override
    public void setParticlePairDistanceParametersCacheActivity(boolean aValue) {
        this.isParticlePairDistanceParametersCacheActive = aValue;
    }

    /**
     * Activity of cache for ParticlePairDistanceParameterss
     * 
     * @return True: Cache for ParticlePairDistanceParameterss is active, false: Otherwise
     */
    @Override
    public boolean getParticlePairDistanceParametersCacheActivity() {
        return this.isParticlePairDistanceParametersCacheActive;
    }
    
    /**
     * Sets cache for ParticlePairDistanceParameterss
     * 
     * @param aParticlePairDistanceParametersCache Cache for ParticlePairDistanceParameterss
     */
    @Override
    public void setParticlePairDistanceParametersCache(ParticlePairDistanceParameters[][] aParticlePairDistanceParametersCache) {
        this.particlePairDistanceParametersCache = aParticlePairDistanceParametersCache;
    }
    
    /**
     * Sets cache for ParticlePairDistanceParameterss
     * 
     * @param aParticlePairInteractionCalculator Particle pair interaction 
     * calculator whose cache is transferred
     */
    @Override
    public void setParticlePairDistanceParametersCache(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) {
        this.particlePairDistanceParametersCache = aParticlePairInteractionCalculator.getParticlePairDistanceParametersCache();
    }

    /**
     * Removes cache for ParticlePairDistanceParameterss
     */
    @Override
    public void removeParticlePairDistanceParametersCache() {
        this.particlePairDistanceParametersCache = null;
    }
    
    /**
     * Returns cache for ParticlePairDistanceParameterss
     * 
     * @return Cache for ParticlePairDistanceParameterss
     */
    @Override
    public ParticlePairDistanceParameters[][] getParticlePairDistanceParametersCache() {
        return this.particlePairDistanceParametersCache;
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Miscellaneous methods">
    /**
     * Executor service shutdown
     */
    @Override
    public void shutdownExecutorService() {
        this.calculatorUtils.shutdownExecutorService();
    }
    
    /**
     * Returns accumulated (total) sum of all potential energy adders
     * (No checks are performed)
     * 
     * @return Accumulated (total) sum of all potential energy adders
     */
    @Override
    public double getAccumulatedPotentialEnergyAddersSum() {
        return this.calculatorUtils.getAccumulatedPotentialEnergyAddersSum();
    }
    
    /**
     * Returns accumulated (total) sum of all pressure tensor diagonal x term adders
     * (No checks are performed)
     * 
     * @return Accumulated (total) sum of all pressure tensor diagonal x term adders
     */
    @Override
    public double getAccumulatedPressureXAddersSum() {
        return this.calculatorUtils.getAccumulatedPressureXAddersSum();
    }
    
    /**
     * Returns accumulated (total) sum of all pressure tensor diagonal y term adders
     * (No checks are performed)
     * 
     * @return Accumulated (total) sum of all pressure tensor diagonal y term adders
     */
    @Override
    public double getAccumulatedPressureYAddersSum() {
        return this.calculatorUtils.getAccumulatedPressureYAddersSum();
    }
    
    /**
     * Returns accumulated (total) sum of all pressure tensor diagonal z term adders
     * (No checks are performed)
     * 
     * @return Accumulated (total) sum of all pressure tensor diagonal z term adders
     */
    @Override
    public double getAccumulatedPressureZAddersSum() {
        return this.calculatorUtils.getAccumulatedPressureZAddersSum();
    }
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Box size
     * 
     * @return Box size
     */
    @Override
    public BoxSize getBoxSize() {
        return this.boxSize;
    }
    
    /**
     * Periodic boundaries
     * 
     * @return Periodic boundaries
     */
    @Override
    public PeriodicBoundaries getPeriodicBoundaries() {
        return this.periodicBoundaries;
    }

    /**
     * Cut-off length for partitioning of the box
     * 
     * @return Cut-off length for partitioning of the box
     */
    @Override
    public double getCutOffLength() {
        return this.cutOffLength;
    }

    /**
     * Parallelisation info
     * 
     * @return Parallelisation info
     */
    @Override
    public ParallelizationInfo getParallelizationInfo() {
        return this.parallelizationInfo;
    }
    
    /**
     * True: Cell-based calculation may be performed, 
     * false: Full loop over particle index pairs is performed
     * 
     * @return True: Cell-based calculation may be performed, 
     * false: Full loop over particle index pairs is performed
     */
    public boolean isCellBasedCalculation() {
        return this.isCellBasedCalculation;
    }

    /**
     * Cell neighbours
     * 
     * @return Cell neighbours or null if isCellBasedCalculation() is false
     */
    @Override
    public int[][] getCellNeigbours() {
        if (this.isCellBasedCalculation) {
            return this.cellNeigbours;
        } else {
            return null;
        }
    }

    /**
     * Cell chunks that are parallelisation-safe
     * 
     * @return Cell chunks that are parallelisation-safe
     */
    @Override
    public int[][] getParallelizationSafeCellChunks() {
        if (this.isCellBasedCalculation) {
            return this.parallelizationSafeCellChunks;
        } else {
            return null;
        }
    }

    /**
     * Head particle index of cell array
     * 
     * @return Head particle index of cell array or null if isCellBasedCalculation() 
     * is false
     */
    @Override
    public int[] getHeadParticleIndexOfCellArray() {
        if (this.isCellBasedCalculation) {
            return this.headParticleIndexOfCell;
        } else {
            return null;
        }
    }
    
    /**
     * Next particle index of cell array
     * 
     * @return Next particle index of cell array or null if isCellBasedCalculation() 
     * is false
     */
    @Override
    public int[] getNextParticleIndexOfCellArray() {
        if (this.isCellBasedCalculation) {
            return this.nextParticleIndexOfCell;
        } else {
            return null;
        }
    }
    
    /**
     * Random number seed
     * 
     * @return Random number seed
     */
    @Override
    public AtomicInteger getRandomNumberSeed() {
        return this.randomNumberSeed;
    }
    
    /**
     * Simulation logger
     * 
     * @return Simulation logger
     */
    @Override
    public ILogger getSimulationLogger() {
        return this.simulationLogger;
    }
    
    /**
     * Factory for new objects
     * 
     * @return Factory for new objects
     */
    @Override
    public Factory getFactory() {
        return this.factory;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (set)">
    /**
     * Set particle assignments (e.g. for call of method
 this.calculateCellBasedParticlePairInteractionsWithoutParticleCellAssignments())
 (No checks are performed)
     * 
     * @param aHeadParticleIndexOfCellArray Head particle index of cell array
     * (retrieved with corresponding getter)
     * @param aNextParticleIndexOfCellArray Next particle index of cell array
     * (retrieved with corresponding getter)
     */
    @Override
    public void setParticleCellAssignments(int[] aHeadParticleIndexOfCellArray, int[] aNextParticleIndexOfCellArray) {
        this.headParticleIndexOfCell = aHeadParticleIndexOfCellArray;
        this.nextParticleIndexOfCell = aNextParticleIndexOfCellArray;
    }

    /**
     * Set particle assignments (e.g. for call of method
 this.calculateCellBasedParticlePairInteractionsWithoutParticleCellAssignments())
 (No checks are performed)
     * 
     * @param aParticlePairInteractionCalculator Particle pair interaction calculator whose 
     * head particle index of cell array and next particle index of cell array are 
     * transferred
     */
    @Override
    public void setParticleCellAssignments(IParticlePairInteractionCalculator aParticlePairInteractionCalculator) {
        this.headParticleIndexOfCell = aParticlePairInteractionCalculator.getHeadParticleIndexOfCellArray();
        this.nextParticleIndexOfCell = aParticlePairInteractionCalculator.getNextParticleIndexOfCellArray();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    // <editor-fold defaultstate="collapsed" desc="- Cell-based methods">
    /**
     * Determines particle index pairs
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    private boolean determineParticleIndexPairs (
        double[] aR_x,
        double[] aR_y,
        double[] aR_z, 
        Parameters aParameters
    ) {
        try {
            int tmpCurrentParallelTaskNumber = this.cellChunkSize/this.parallelizationInfo.getMinimumParallelTaskCellNumber() + 1;
            if (tmpCurrentParallelTaskNumber > this.parallelizationInfo.getParallelTaskNumber()) {
                tmpCurrentParallelTaskNumber = this.parallelizationInfo.getParallelTaskNumber();
            }
            this.parallelizationInfo.setMaximumUsedParallelTaskNumber(tmpCurrentParallelTaskNumber);
            if (this.isParticlePairDistanceParametersCacheActive) {
                if (this.particlePairDistanceParametersCache == null) {
                    this.particlePairDistanceParametersCache = new ParticlePairDistanceParameters[this.parallelizationSafeCellChunks.length][tmpCurrentParallelTaskNumber];
                    int tmpCapacity = aR_x.length/(this.parallelizationSafeCellChunks.length * tmpCurrentParallelTaskNumber);
                    for (ParticlePairDistanceParameters[] tmpParticlePairDistanceParametersCachePerCellChunk : this.particlePairDistanceParametersCache) {
                        for (int i = 0; i < tmpParticlePairDistanceParametersCachePerCellChunk.length; i++) {
                            tmpParticlePairDistanceParametersCachePerCellChunk[i] = new ParticlePairDistanceParameters(tmpCapacity);
                        }
                    }
                } else {
                    this.resetParticlePairDistanceParametersCache();
                }
            }
            if (this.parallelizationInfo.getParallelTaskNumber() > 1 && tmpCurrentParallelTaskNumber > 1) {
                // <editor-fold defaultstate="collapsed" desc="Parallelized calculation">
                // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                this.simulationLogger.appendParallelization("ParticlePairInteractionCalculator.determineParticleIndexPairs: Current parallel task number = " + String.valueOf(tmpCurrentParallelTaskNumber));
                // </editor-fold>
                if (this.calculatorUtils.getAdderGroups() == null) {
                    RandomAdderGroup[] tmpRandomAdderGroups = new RandomAdderGroup[tmpCurrentParallelTaskNumber];
                    for (int i = 0; i < tmpCurrentParallelTaskNumber; i++) {
                        tmpRandomAdderGroups[i] = new RandomAdderGroup(this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet()));
                    }
                    this.calculatorUtils.setAdderGroups(tmpRandomAdderGroups);
                }
                if (this.calculatorUtils.getExecutorService() == null) {
                    this.calculatorUtils.setExecutorService(Executors.newFixedThreadPool(tmpCurrentParallelTaskNumber));
                }
                RandomAdderGroup[] tmpRandomAdderGroups = ((RandomAdderGroup[]) this.calculatorUtils.getAdderGroups());
                ParticlePairDistanceParameters[] tmpParticlePairDistanceParametersCachePerCellChunk = null;
                ParticlePairDistanceParameters tmpParticlePairDistanceParameters = null;
                for (int i = 0; i < this.parallelizationSafeCellChunks.length; i++) {
                    int tmpCellNumberForParallelizedTask = this.cellChunkSize/tmpCurrentParallelTaskNumber;
                    int[] tmpParallelizationSafeCellChunk = this.parallelizationSafeCellChunks[i];
                    LinkedList<CellChunkCalculationTask> tmpTaskList = new LinkedList<>();
                    int tmpStartIndex = 0;
                    int tmpEndIndex = tmpCellNumberForParallelizedTask;
                    if (this.isParticlePairDistanceParametersCacheActive) {
                        tmpParticlePairDistanceParametersCachePerCellChunk = this.particlePairDistanceParametersCache[i];
                    }
                    // Loop exclusive the last task ...
                    for (int k = 0; k < tmpCurrentParallelTaskNumber - 1; k++) {
                        if (this.isParticlePairDistanceParametersCacheActive) {
                            tmpParticlePairDistanceParameters = tmpParticlePairDistanceParametersCachePerCellChunk[k];
                        }
                        CellChunkCalculationTask newTask = 
                            new CellChunkCalculationTask(
                                this.simulationLogger,
                                tmpParallelizationSafeCellChunk, 
                                tmpStartIndex, 
                                tmpEndIndex, 
                                aR_x, 
                                aR_y, 
                                aR_z, 
                                tmpRandomAdderGroups[k],
                                tmpParticlePairDistanceParameters,
                                aParameters,
                                this);
                        tmpTaskList.add(newTask);
                        tmpStartIndex = tmpEndIndex;
                        tmpEndIndex += tmpCellNumberForParallelizedTask;
                    }
                    // ... to get a correct final end index:
                    if (this.isParticlePairDistanceParametersCacheActive) {
                        tmpParticlePairDistanceParameters = this.particlePairDistanceParametersCache[i][tmpCurrentParallelTaskNumber - 1];
                    }
                    tmpTaskList.add(
                        new CellChunkCalculationTask(
                            this.simulationLogger,
                            tmpParallelizationSafeCellChunk, 
                            tmpStartIndex, 
                            tmpParallelizationSafeCellChunk.length, 
                            aR_x, 
                            aR_y, 
                            aR_z, 
                            tmpRandomAdderGroups[tmpCurrentParallelTaskNumber - 1],
                            tmpParticlePairDistanceParameters,
                            aParameters, 
                            this
                        )
                    );
                    this.calculatorUtils.getExecutorService().invokeAll(tmpTaskList);
                }
                // </editor-fold>
            } else {
                // <editor-fold defaultstate="collapsed" desc="Sequential calculation">
                // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                this.simulationLogger.appendParallelization("ParticlePairInteractionCalculator.determineParticleIndexPairs: Sequential calculation.");
                // </editor-fold>
                if (this.calculatorUtils.getAdderGroups() == null) {
                    this.calculatorUtils.setAdderGroups(new RandomAdderGroup[] {new RandomAdderGroup(this.factory.getNewOrJumpedRandomNumberGenerator(this.randomNumberSeed.incrementAndGet()))});
                }
                RandomAdderGroup tmpRandomAdderGroup = ((RandomAdderGroup[]) this.calculatorUtils.getAdderGroups())[0];
                ParticlePairDistanceParameters tmpParticlePairDistanceParameters = null;
                for (int i = 0; i < this.parallelizationSafeCellChunks.length; i++) {
                    int[] tmpParallelizationSafeCellChunk = this.parallelizationSafeCellChunks[i];
                    if (this.isParticlePairDistanceParametersCacheActive) {
                        tmpParticlePairDistanceParameters = this.particlePairDistanceParametersCache[i][0];
                    }
                    for (int tmpCellIndex = 0; tmpCellIndex < tmpParallelizationSafeCellChunk.length; tmpCellIndex++) {
                        this.calculateSingleCellParticlePairInteractions(
                            tmpParallelizationSafeCellChunk[tmpCellIndex],
                            aR_x,
                            aR_y,
                            aR_z,
                            tmpRandomAdderGroup,
                            tmpParticlePairDistanceParameters,
                            aParameters
                        );
                    }
                }
                // </editor-fold>
            }
            return true;
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("ParticlePairInteractionCalculator.determineParticleIndexPairs", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
    }
    /**
     * Support method
     * (No checks are performed)
     * 
     * @param aCellIndex Cell index
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aRandomAdderGroup Random adder group
     * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance for cache (may be null)
     * @param aParameters Parameters (may be null)
     */
    private void calculateSingleCellParticlePairInteractions(
        int aCellIndex,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        RandomAdderGroup aRandomAdderGroup,
        ParticlePairDistanceParameters aParticlePairDistanceParameters,
        Parameters aParameters
    ) {
        int tmpParticleIndex1 = this.headParticleIndexOfCell[aCellIndex];
        while (tmpParticleIndex1 != -1) {
            // <editor-fold defaultstate="collapsed" desc="Particle index pairs with particle indices of same cell">
            int tmpParticleIndex2 = this.nextParticleIndexOfCell[tmpParticleIndex1];
            while (tmpParticleIndex2 != -1) {
                this.evaluateParticlePairInteraction(
                    tmpParticleIndex1, 
                    tmpParticleIndex2, 
                    aR_x, 
                    aR_y, 
                    aR_z, 
                    aRandomAdderGroup,
                    aParticlePairDistanceParameters,
                    aParameters
                );
                tmpParticleIndex2 = this.nextParticleIndexOfCell[tmpParticleIndex2];
            }
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Particle index pairs with particle indices of neighbour cells">
            int[] tmpNeighbours = this.cellNeigbours[aCellIndex];
            for (int tmpNeighbourIndex = 0; tmpNeighbourIndex < tmpNeighbours.length; tmpNeighbourIndex++) {
                int tmpNeighbourCellIndex = tmpNeighbours[tmpNeighbourIndex];
                if (tmpNeighbourCellIndex == -1) {
                    // Edge cells may not necessarily have all neighbours due to possible non-periodic boundaries
                    break;
                } else {
                    tmpParticleIndex2 = this.headParticleIndexOfCell[tmpNeighbourCellIndex];
                    while (tmpParticleIndex2 != -1) {
                        this.evaluateParticlePairInteraction(
                            tmpParticleIndex1, 
                            tmpParticleIndex2, 
                            aR_x, 
                            aR_y, 
                            aR_z, 
                            aRandomAdderGroup,
                            aParticlePairDistanceParameters,
                            aParameters
                        );
                        tmpParticleIndex2 = this.nextParticleIndexOfCell[tmpParticleIndex2];
                    }
                }
            }
            // </editor-fold>
            tmpParticleIndex1 = this.nextParticleIndexOfCell[tmpParticleIndex1];
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Miscellaneous methods">
    /**
     * Returns index along x-direction
     * 
     * @param aPosition Position along x-direction
     * @return Index along x-direction
     */
    private int getIndexAlongX(double aPosition) {
        // Fast implementation for positive values:
        return (int) (aPosition/this.cutOffLengthX);
        // Slower alternative:
        // return (int) FastMath.floor(aPosition/this.cutOffLengthX);
    }

    /**
     * Returns index along y-direction
     * 
     * @param aPosition Position along y-direction
     * @return Index along y-direction
     */
    private int getIndexAlongY(double aPosition) {
        // Fast implementation for positive values:
        return (int) (aPosition/this.cutOffLengthY);
        // Slower alternative:
        // return (int) FastMath.floor(aPosition/this.cutOffLengthY);
    }

    /**
     * Returns index along z-direction
     * 
     * @param aPosition Position along z-direction
     * @return Index along z-direction
     */
    private int getIndexAlongZ(double aPosition) {
        // Fast implementation for positive values:
        return (int) (aPosition/this.cutOffLengthZ);
        // Slower alternative:
        // return (int) FastMath.floor(aPosition/this.cutOffLengthZ);
    }
    
    /**
     * Support method
     * (No checks are performed)
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aRandomAdderGroup Random adder group (may be null)
     * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance for cache (may be null)
     * @param aParameters Parameters (may be null)
     */
    private void evaluateParticlePairInteraction(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        RandomAdderGroup aRandomAdderGroup,
        ParticlePairDistanceParameters aParticlePairDistanceParameters,
        Parameters aParameters
    ) {
        double tmpCorrectedRij_x = 
            Utils.correctPositionDifference(
                aR_x[aParticleIndex_i] - aR_x[aParticleIndex_j],
                this.periodicBoundaries.isAlongX(),
                this.boxSize.getXHalfLength(),
                this.boxSize.getXLength(),
                this.boxSize.getNegativeXHalfLength()
            );
        double tmpCorrectedRij_x_Square = tmpCorrectedRij_x * tmpCorrectedRij_x;
        if (tmpCorrectedRij_x_Square < this.cutOffLengthSquare) {
            double tmpCorrectedRij_y = 
                Utils.correctPositionDifference(
                    aR_y[aParticleIndex_i] - aR_y[aParticleIndex_j],
                    this.periodicBoundaries.isAlongY(),
                    this.boxSize.getYHalfLength(),
                    this.boxSize.getYLength(),
                    this.boxSize.getNegativeYHalfLength()
                );
            double tmpCorrectedRij_y_Square = tmpCorrectedRij_y * tmpCorrectedRij_y;
            if (tmpCorrectedRij_y_Square < this.cutOffLengthSquare) {
                double tmpCorrectedRij_z = 
                    Utils.correctPositionDifference(
                        aR_z[aParticleIndex_i] - aR_z[aParticleIndex_j],
                        this.periodicBoundaries.isAlongZ(),
                        this.boxSize.getZHalfLength(),
                        this.boxSize.getZLength(),
                        this.boxSize.getNegativeZHalfLength()
                    );
                double tmpCorrectedRij_z_Square = tmpCorrectedRij_z * tmpCorrectedRij_z;
                if (tmpCorrectedRij_z_Square < this.cutOffLengthSquare) {
                    double tmpRij_Square = tmpCorrectedRij_x_Square + tmpCorrectedRij_y_Square + tmpCorrectedRij_z_Square;
                    if (tmpRij_Square < this.cutOffLengthSquare) {
                        this.calculateParticlePairInteraction(
                            aParticleIndex_i, 
                            aParticleIndex_j, 
                            tmpCorrectedRij_x,
                            tmpCorrectedRij_y,
                            tmpCorrectedRij_z,
                            tmpRij_Square, 
                            aRandomAdderGroup,
                            aParameters,
                            aParticlePairDistanceParameters
                        );
                    }
                }
            }
        }
    }
    // <editor-fold defaultstate="collapsed" desc="- Cache related methods">
    /**
     * Resets this.particlePairDistanceParametersCache
     * NOTE: No checks are performed!
     */
    private void resetParticlePairDistanceParametersCache() {
        for (ParticlePairDistanceParameters[] tmpParticlePairDistanceParametersCachePerCellChunk : this.particlePairDistanceParametersCache) {
            for (ParticlePairDistanceParameters tmpParticlePairDistanceParametersCachePerParallelTask : tmpParticlePairDistanceParametersCachePerCellChunk) {
                tmpParticlePairDistanceParametersCachePerParallelTask.reset();
            }
        }
    }
    // </editor-fold>
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Abstract protected methods">
    /**
     * Abstract method for particle pair interaction calculation 
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or double adder 
     * is necessary.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aRij_Square Squared distance between particle i and j
     * @param aRandomAdderGroup Random adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     * @param aParticlePairDistanceParameters ParticlePairDistanceParameters instance ((may be null)
     */
    protected abstract void calculateParticlePairInteraction(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters,
        ParticlePairDistanceParameters aParticlePairDistanceParameters
    );

    /**
     * Abstract method for particle pair interaction calculation 
     * NOTE: ParticlePairInteractionCalculator parallelisation guarantees that
     * NO thread-safe implementation of random number generator or double adder 
     * is necessary.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aRij_Square Squared distance between particle i and j
     * @param aRij Distance between particle i and j
     * @param aRandomAdderGroup Random adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     */
    protected abstract void calculateParticlePairInteraction(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aRij_Square,
        double aRij,
        RandomAdderGroup aRandomAdderGroup,
        Parameters aParameters);
    // </editor-fold>

}
