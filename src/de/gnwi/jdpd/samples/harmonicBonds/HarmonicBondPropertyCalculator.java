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
package de.gnwi.jdpd.samples.harmonicBonds;

import de.gnwi.jdpd.interfaces.IHarmonicBondPropertyCalculator;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.utilities.CalculatorUtils;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.utilities.AdderGroup;
import de.gnwi.jdpd.utilities.Utils;
import java.util.LinkedList;
import java.util.concurrent.Callable;
import java.util.concurrent.Executors;

/**
 * Harmonic bond property calculator for parallelised calculation of bond chunk 
 * properties.
 * 
 * The bond chunk parallelisation implementation allows the use of fast NON 
 * thread-safe objects and accumulation items like non thread-safe double adders 
 * and standard (NON atomic) double arrays.
 * 
 * @author Achim Zielesny
 */
public abstract class HarmonicBondPropertyCalculator implements IHarmonicBondPropertyCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private class PartialBondChunkCalculationTask">
    /**
     * Callable for calculation of partial bond chunks in a parallelised manner 
     * with executor service
     */
    private class PartialBondChunkCalculationTask implements Callable<Boolean> {
        
        // <editor-fold defaultstate="collapsed" desc="Private final class variables">
        /**
         * Bond chunk arrays that are parallelisation-safe
         */
        private final HarmonicBondChunkArrays parallelizationSafeBondChunkArrays;

        /**
         * Start index for parallelizationSafeBondChunkArrays
         */
        private final int startIndex;

        /**
         * End index for parallelizationSafeBondChunkArrays
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
         * Adder group
         */
        private final AdderGroup adderGroup;
        
        /**
         * Particle pair interaction calculator
         */
        private final HarmonicBondPropertyCalculator bondPropertyCalculator;
        // </editor-fold>
        //
        // <editor-fold defaultstate="collapsed" desc="Constructor">
        /**
         * Constructor
         * 
         * @param aParallelizationSafeBondChunkArrays Bond chunk arrays that are parallelisation-safe
         * @param aStartIndex Start index for parallelizationSafeCellChunk
         * @param anEndIndex End index for parallelizationSafeCellChunk
         * @param aR_x Current x-components of particle positions in simulation box
         * @param aR_y Current y-components of particle positions in simulation box
         * @param aR_z Current z-components of particle positions in simulation box
         * @param anAdderGroup Adder group (NOTE: NOT thread-safe)
         * @param aParameters Parameters (may be null)
         * @param aHarmonicBondPropertyCalculator Bond property calculator
         */
        public PartialBondChunkCalculationTask(
            HarmonicBondChunkArrays aParallelizationSafeBondChunkArrays,
            int aStartIndex, 
            int anEndIndex,
            double[] aR_x,
            double[] aR_y,
            double[] aR_z,
            AdderGroup anAdderGroup,
            Parameters aParameters,
            HarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) {
            this.parallelizationSafeBondChunkArrays = aParallelizationSafeBondChunkArrays;
            this.startIndex = aStartIndex;
            this.endIndex = anEndIndex;
            this.r_x = aR_x;
            this.r_y = aR_y;
            this.r_z = aR_z;
            this.adderGroup = anAdderGroup;
            this.parameters = aParameters;
            this.bondPropertyCalculator = aHarmonicBondPropertyCalculator;
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
            int[] tmpParticleIndices1 = this.parallelizationSafeBondChunkArrays.getParticleIndices1();
            int[] tmpParticleIndices2 = this.parallelizationSafeBondChunkArrays.getParticleIndices2();
            double[] tmpBondLengths = this.parallelizationSafeBondChunkArrays.getBondLengths();
            double[] tmpForceConstants = this.parallelizationSafeBondChunkArrays.getForceConstants();
            boolean[] tmpRepulsionFlags = this.parallelizationSafeBondChunkArrays.getRepulsionFlags();
            for (int i = this.startIndex; i < this.endIndex; i++) {
                this.bondPropertyCalculator.correctAndCalculateBondProperty(
                    tmpParticleIndices1[i],
                    tmpParticleIndices2[i],
                    tmpBondLengths[i],
                    tmpForceConstants[i],
                    tmpRepulsionFlags[i],
                    this.r_x,
                    this.r_y,
                    this.r_z,
                    this.adderGroup,
                    this.parameters
                );
            }
            return true;
        }
        // </editor-fold>
        
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Calculator utils
     */
    private final CalculatorUtils calculatorUtils;
    
    /**
     * Box size
     */
    private final BoxSize boxSize;
    
    /**
     * Periodic boundaries
     */
    private final PeriodicBoundaries periodicBoundaries;
    
    /**
     * Parallelisation info
     */
    private final ParallelizationInfo parallelizationInfo;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected final class variables">
    /**
     * Simulation simulationLogger
     */
    protected final ILogger simulationLogger;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation simulationLogger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aParallelizationInfo Parallelisation info
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected HarmonicBondPropertyCalculator(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        ParallelizationInfo aParallelizationInfo) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("HarmonicBondPropertyCalculator.Constructor: aSimulationLogger is null.");
        }
        if (aBoxSize == null) {
            throw new IllegalArgumentException("HarmonicBondPropertyCalculator.Constructor: aBoxSize is null.");
        }
        if (aPeriodicBoundaries == null) {
            throw new IllegalArgumentException("HarmonicBondPropertyCalculator.Constructor: aPeriodicBoundaries is null.");
        }
        if (aParallelizationInfo == null) {
            throw new IllegalArgumentException("HarmonicBondPropertyCalculator.Constructor: aParallelizationInfo is null.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        try {
            this.boxSize = aBoxSize;
            this.periodicBoundaries = aPeriodicBoundaries;
            this.parallelizationInfo = aParallelizationInfo;
            this.calculatorUtils = new CalculatorUtils();
            // <editor-fold defaultstate="collapsed" desc="Method call logging">
            this.simulationLogger.appendMethodCall("HarmonicBondPropertyCalculator.Constructor");
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("HarmonicBondPropertyCalculator.Constructor", Utils.getStacktrace(anException));
            // </editor-fold>
            throw anException;
        }
    }
    
    /**
     * Constructor that clones aHarmonicBondPropertyCalculator
     * 
     * @param aHarmonicBondPropertyCalculator BondPropertyCalculator instance
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected HarmonicBondPropertyCalculator(IHarmonicBondPropertyCalculator aHarmonicBondPropertyCalculator) throws IllegalArgumentException {
        this(
            aHarmonicBondPropertyCalculator.getSimulationLogger(), 
            aHarmonicBondPropertyCalculator.getBoxSize(), 
            aHarmonicBondPropertyCalculator.getPeriodicBoundaries(), 
            aHarmonicBondPropertyCalculator.getParallelizationInfo()
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("HarmonicBondPropertyCalculator.Constructor WITH aHarmonicBondPropertyCalculator");
        // </editor-fold>
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    // <editor-fold defaultstate="collapsed" desc="- Calculation methods">
    /**
     * Calculates bond properties.
     * NOTE: NO checks are performed.
     * 
     * @param aBondChunkArraysList List of bond chunk arrays
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aParameters Parameters (may be null)
     * @return True: Operation successful, false: Otherwise
     */
    @Override
    public boolean calculateBondProperties(
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        Parameters aParameters) {
        try {
            // <editor-fold defaultstate="collapsed" desc="Reset double adders">
            this.calculatorUtils.resetAdders();
            // </editor-fold>
            // <editor-fold defaultstate="collapsed" desc="Calculate bond properties">
            for (HarmonicBondChunkArrays tmpBondChunkArrays : aBondChunkArraysList) {
                int tmpCurrentParallelTaskNumber = tmpBondChunkArrays.getLength()/this.parallelizationInfo.getMinimumParallelTaskHarmonicBondNumber() + 1;
                if (this.parallelizationInfo.getParallelTaskNumber() > 1 && tmpCurrentParallelTaskNumber > 1) {
                    // <editor-fold defaultstate="collapsed" desc="Parallelized calculation">
                    if (tmpCurrentParallelTaskNumber > this.parallelizationInfo.getParallelTaskNumber()) {
                        tmpCurrentParallelTaskNumber = this.parallelizationInfo.getParallelTaskNumber();
                    }
                    this.parallelizationInfo.setMaximumUsedParallelTaskNumber(tmpCurrentParallelTaskNumber);
                    // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                    this.simulationLogger.appendParallelization("HarmonicBondPropertyCalculator.calculateBondProperties: Current parallel task number = " + String.valueOf(tmpCurrentParallelTaskNumber));
                    // </editor-fold>
                    if (this.calculatorUtils.getAdderGroups() == null) {
                        AdderGroup[] tmpAdderGroups = new AdderGroup[tmpCurrentParallelTaskNumber];
                        for (int i = 0; i < tmpCurrentParallelTaskNumber; i++) {
                            tmpAdderGroups[i] = new AdderGroup();
                        }
                        this.calculatorUtils.setAdderGroups(tmpAdderGroups);
                    }
                    AdderGroup[] tmpAdderGroups = this.calculatorUtils.getAdderGroups();
                    int tmpBondNumberForParallelizedTask = tmpBondChunkArrays.getLength()/tmpCurrentParallelTaskNumber;
                    LinkedList<PartialBondChunkCalculationTask> tmpTaskList = new LinkedList<>();
                    int tmpStartIndex = 0;
                    int tmpEndIndex = tmpBondNumberForParallelizedTask;
                    // Loop exclusive the last task ...
                    for (int k = 0; k < tmpCurrentParallelTaskNumber - 1; k++) {
                        PartialBondChunkCalculationTask newTask = 
                            new PartialBondChunkCalculationTask(
                                tmpBondChunkArrays,
                                tmpStartIndex, 
                                tmpEndIndex, 
                                aR_x,
                                aR_y,
                                aR_z,
                                tmpAdderGroups[k],
                                aParameters,    
                                this);
                        tmpTaskList.add(newTask);
                        tmpStartIndex = tmpEndIndex;
                        tmpEndIndex += tmpBondNumberForParallelizedTask;
                    }
                    // ... to get a correct final end index:
                    tmpTaskList.add(
                        new PartialBondChunkCalculationTask(
                            tmpBondChunkArrays, 
                            tmpStartIndex, 
                            tmpBondChunkArrays.getLength(), 
                            aR_x,
                            aR_y,
                            aR_z,
                            tmpAdderGroups[tmpCurrentParallelTaskNumber - 1],
                            aParameters,    
                            this
                        )
                    );
                    if (this.calculatorUtils.getExecutorService() == null) {
                        this.calculatorUtils.setExecutorService(Executors.newFixedThreadPool(tmpCurrentParallelTaskNumber));
                    }
                    this.calculatorUtils.getExecutorService().invokeAll(tmpTaskList);
                    // </editor-fold>
                } else {
                    // <editor-fold defaultstate="collapsed" desc="Sequential calculation">
                    // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                    this.simulationLogger.appendParallelization("HarmonicBondPropertyCalculator.calculateBondProperties: Sequential calculation.");
                    // </editor-fold>
                    if (this.calculatorUtils.getAdderGroups() == null) {
                        this.calculatorUtils.setAdderGroups(new AdderGroup[] {new AdderGroup()});
                    }
                    AdderGroup tmpAdderGroup = this.calculatorUtils.getAdderGroups()[0];
                    int[] tmpParticleIndices1 = tmpBondChunkArrays.getParticleIndices1();
                    int[] tmpParticleIndices2 = tmpBondChunkArrays.getParticleIndices2();
                    double[] tmpBondLengths = tmpBondChunkArrays.getBondLengths();
                    double[] tmpForceConstants = tmpBondChunkArrays.getForceConstants();
                    boolean[] tmpRepulsionFlags = tmpBondChunkArrays.getRepulsionFlags();
                    for (int i = 0; i < tmpParticleIndices1.length; i++) {
                        this.correctAndCalculateBondProperty(
                            tmpParticleIndices1[i],
                            tmpParticleIndices2[i],
                            tmpBondLengths[i],
                            tmpForceConstants[i],
                            tmpRepulsionFlags[i],
                            aR_x,
                            aR_y,
                            aR_z,
                            tmpAdderGroup,
                            aParameters
                        );
                    }
                    // </editor-fold>
                }
            }
            return true;
            // </editor-fold>
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("HarmonicBondPropertyCalculator.calculateBondProperties", Utils.getStacktrace(anException));
            // </editor-fold>
            return false;
        }
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
     * Return accumulated (total) sum of all potential energy adders
     * 
     * @return Return accumulated (total) sum of all potential energy adders
     */
    @Override
    public double getAccumulatedPotentialEnergyAddersSum() {
        return this.calculatorUtils.getAccumulatedPotentialEnergyAddersSum();
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal x term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal x term adders
     */
    @Override
    public double getAccumulatedPressureXAddersSum() {
        return this.calculatorUtils.getAccumulatedPressureXAddersSum();
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal y term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal y term adders
     */
    @Override
    public double getAccumulatedPressureYAddersSum() {
        return this.calculatorUtils.getAccumulatedPressureYAddersSum();
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal z term adders
     * NOTE: NO checks are performed.
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal z term adders
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
     * Parallelisation info
     * 
     * @return Parallelisation info
     */
    @Override
    public ParallelizationInfo getParallelizationInfo() {
        return this.parallelizationInfo;
    }
    
    /**
     * Simulation simulationLogger
     * 
     * @return Simulation simulationLogger
     */
    @Override
    public ILogger getSimulationLogger() {
        return this.simulationLogger;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Support method
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aBondLength Bond length
     * @param aForceConstant Force constant
     * @param anIsRepulsion True: Repulsion for bond is to be calculated, false: Otherwise (no repulsion, attraction only)
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param anAdderGroup Adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     */
    private void correctAndCalculateBondProperty(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aBondLength,
        double aForceConstant,
        boolean anIsRepulsion,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        AdderGroup anAdderGroup,
        Parameters aParameters) {
        double tmpCorrectedRij_x = 
            Utils.correctPositionDifference(
                aR_x[aParticleIndex_i] - aR_x[aParticleIndex_j],
                this.periodicBoundaries.isAlongX(),
                this.boxSize.getXHalfLength(),
                this.boxSize.getXLength(),
                this.boxSize.getNegativeXHalfLength()
            );
        double tmpCorrectedRij_y = 
            Utils.correctPositionDifference(
                aR_y[aParticleIndex_i] - aR_y[aParticleIndex_j],
                this.periodicBoundaries.isAlongY(),
                this.boxSize.getYHalfLength(),
                this.boxSize.getYLength(),
                this.boxSize.getNegativeYHalfLength()
            );
        double tmpCorrectedRij_z = 
            Utils.correctPositionDifference(
                aR_z[aParticleIndex_i] - aR_z[aParticleIndex_j],
                this.periodicBoundaries.isAlongZ(),
                this.boxSize.getZHalfLength(),
                this.boxSize.getZLength(),
                this.boxSize.getNegativeZHalfLength()
            );
        this.calculateBondProperty(
            aParticleIndex_i, 
            aParticleIndex_j, 
            tmpCorrectedRij_x,
            tmpCorrectedRij_y,
            tmpCorrectedRij_z,
            aBondLength, 
            aForceConstant,
            anIsRepulsion,
            anAdderGroup,
            aParameters
        );
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Abstract protected methods">
    /**
     * Abstract method for bond property calculation
     * NOTE: HarmonicBondPropertyCalculator parallelisation guarantees for 
     * aDoubleAdder being unique in current working thread thus NO thread-safe 
     * implementation of double adder is necessary.
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x x[aParticleIndex_i] - x[aParticleIndex_j] 
     * @param aRij_y y[aParticleIndex_i] - y[aParticleIndex_j] 
     * @param aRij_z z[aParticleIndex_i] - z[aParticleIndex_j] 
     * @param aBondLength Bond length
     * @param aForceConstant Force constant
     * @param anIsRepulsion True: Repulsion for bond is to be calculated, false: Otherwise (no repulsion, attraction only)
     * @param anAdderGroup Adder group (NOTE: NOT thread-safe)
     * @param aParameters Parameters (may be null)
     */
    protected abstract void calculateBondProperty(
        int aParticleIndex_i, 
        int aParticleIndex_j, 
        double aRij_x,
        double aRij_y,
        double aRij_z,
        double aBondLength,
        double aForceConstant,
        boolean anIsRepulsion,
        AdderGroup anAdderGroup,
        Parameters aParameters
    );
    // </editor-fold>

}
