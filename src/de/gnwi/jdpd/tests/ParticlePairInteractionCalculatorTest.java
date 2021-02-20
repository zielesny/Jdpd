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
package de.gnwi.jdpd.tests;

import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.samples.logger.MemoryLogger;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpd.parameters.ParticleArrays;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.parameters.Parameters;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import de.gnwi.jdpd.parameters.TestObjects;
import de.gnwi.jdpd.samples.random.ApacheCommonsRandom;
import de.gnwi.jdpd.tests.interactions.ParticleIndexPairCalculatorForTesting;
import de.gnwi.jdpd.tests.interactions.ParticleIndexPairNumberCalculatorForTesting;
import de.gnwi.jdpd.tests.interactions.ParticleIndexPairRandomCalculatorForTesting;
import de.gnwi.jdpd.utilities.ParticlePairDistanceParameters;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import junit.framework.*;
import org.apache.commons.rng.simple.RandomSource;


/**
 * Test class for ParticlePairInteractionCalculatorForTesting
 * 
 * @author Achim Zielesny
 */
public class ParticlePairInteractionCalculatorTest extends TestCase {

    // <editor-fold defaultstate="collapsed" desc="Public tests">
    /**
     * Test method
     */
    public void test_IndexPairGeneration_Sequential() {
        this.test_IndexPairGeneration(1);
    }

    /**
     * Test method
     */
    public void test_IndexPairGeneration_Parallel() {
        this.test_IndexPairGeneration(8);
    }

    /**
     * Test method
     */
    public void test_ParticlePairDistanceParametersCache_Sequential() {
        this.test_ParticlePairDistanceParametersCache(1);
    }
    
    /**
     * Test method
     */
    public void test_ParticlePairDistanceParametersCache_Parallel() {
        this.test_ParticlePairDistanceParametersCache(8);
    }

    /**
     * Test method
     */
    public void test_IndexPairRandomGeneration_Sequential() {
        this.test_IndexPairRandomGeneration(1);
    }
    
    /**
     * Test method
     */
    public void test_IndexPairRandomGeneration_Parallel() {
        this.test_IndexPairRandomGeneration(8);
    }

    /**
     * Test method
     */
    public void test_IndexPairCount_Sequential() {
        this.test_IndexPairCount(1);
    }
    
    /**
     * Test method
     */
    public void test_IndexPairCount_Parallel() {
        this.test_IndexPairCount(8);
    }

    /**
     * Test method
     */
    public void test_IndexPairCount_ParallelSize_Performance() {
        long tmpStart;
        long tmpEnd;
        long tmpNumber1;

        System.out.println("test_IndexPairCount_ParallelSize_Performance");
        
        int tmpParallelTaskNumber = 16;
        String tmpOffset = "> ";
        for (double tmpBoxLength = 20.0; tmpBoxLength <= 200.0; tmpBoxLength += 10.0) {
            // NOTE for tests:
            // tmpCutOffLength = 1: Box 20 x 20 x 20 = 8.000, Density 3 * 8000 = 24.000 particles
            // tmpCutOffLength = 1: Box 100 x 100 x 100 = 1.000.000, Density 3 * 1.000.000 = 3.000.000 particles
            System.out.println("Box length      = " + String.valueOf(tmpBoxLength));
            double tmpCutOffLength = 1.0;
            // double tmpBoxLength = 100.0;
            int tmpParticleNumber = (int) (tmpBoxLength * tmpBoxLength * tmpBoxLength * 3);
            System.out.println("Particle number = " + String.valueOf(tmpParticleNumber));
            AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
            BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
            PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

            MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
            ParticleArrays tmpParticles = new ParticleArrays(
                new double[tmpParticleNumber],
                new double[tmpParticleNumber],
                new double[tmpParticleNumber],
                null, null, null, null, 
                null, null, null, null, null, null, null, null,
                null, null, null, null);
            ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
            for (int i = 0; i < tmpParticleNumber; i++) {
                tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
                tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
                tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            }

            // new ParticleIndexPairNumberCalculator -------------------------------

            tmpStart = System.currentTimeMillis();
            ParticleIndexPairNumberCalculatorForTesting tmpParticleIndexPairNumberCalculator = new ParticleIndexPairNumberCalculatorForTesting(new Factory(Factory.RandomType.ACRNG_MT_64, 1, Factory.DpdType.CUTOFF_LENGTH_ONE, Factory.ElectrostaticsType.AD_HOC, Factory.BondType.HARMONIC, Factory.IntegrationType.GWMVV, new Double[] {0.65}), tmpMemoryLogger, tmpBoxSize, tmpPeriodicBoundaries, tmpCutOffLength, null, null, new ParallelizationInfo(2, 2, tmpParallelTaskNumber), tmpRandomNumberSeed);
            tmpEnd = System.currentTimeMillis();
            System.out.println(tmpOffset + "new ParticleIndexPairNumberCalculator,                  t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));

            // Parallel cell-based ---------------------------------------------

            TestObjects tmpTestObjects = new TestObjects(
                new ConcurrentLinkedQueue<>(),
                new ConcurrentLinkedQueue<>(),
                new AtomicInteger(0),
                null
            );
            Parameters tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
            tmpStart = System.currentTimeMillis();
            tmpParticleIndexPairNumberCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
                tmpParticles.getR_x(), 
                tmpParticles.getR_y(), 
                tmpParticles.getR_z(), 
                tmpParameters
            );
            tmpEnd = System.currentTimeMillis();
            System.out.println(tmpOffset + "calculateCellBasedParticlePairInteractions(parallel),   t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
            tmpNumber1 = tmpParameters.getTestObjects().getParticleIndexPairCounter().intValue();
            System.out.println(tmpOffset + "Number of particle index pairs                                = " + String.valueOf(tmpNumber1));
            System.out.println(tmpOffset + "Time in nanoseconds per particle index pair,           t [ns] = " + String.valueOf((int) ((double) (tmpEnd - tmpStart)/(double) tmpNumber1 * 1000000.0)));

            assertTrue("Parallel cell-based", tmpNumber1 > 0);

            System.out.println("");
        }        
    }

    /**
     * Test method
     */
    public void test_IndexPairCount_ParallelProcessor_Performance() {
        long tmpStart;
        long tmpEnd;
        long tmpNumber1;

        System.out.println("test_IndexPairCount_ParallelProcessor_Performance");
        
        String tmpOffset = "> ";
        for (int tmpParallelTaskNumber = 16; tmpParallelTaskNumber >= 1; tmpParallelTaskNumber--) {
            // NOTE for tests:
            // tmpCutOffLength = 1: Box 20 x 20 x 20 = 8.000, Density 3 * 8000 = 24.000 particles
            // tmpCutOffLength = 1: Box 100 x 100 x 100 = 1.000.000, Density 3 * 1.000.000 = 3.000.000 particles
            System.out.println("Number of processors = " + String.valueOf(tmpParallelTaskNumber));
            int tmpBoxLength = 150;
            System.out.println("Box length           = " + String.valueOf(tmpBoxLength));
            double tmpCutOffLength = 1.0;
            // double tmpBoxLength = 100.0;
            int tmpParticleNumber = (int) (tmpBoxLength * tmpBoxLength * tmpBoxLength * 3);
            System.out.println("Particle number      = " + String.valueOf(tmpParticleNumber));
            AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
            BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
            PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

            MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
            ParticleArrays tmpParticles = new ParticleArrays(
                new double[tmpParticleNumber],
                new double[tmpParticleNumber],
                new double[tmpParticleNumber],
                null, null, null, null,
                null, null, null, null, null, null, null, null,
                null, null, null, null);
            ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
            for (int i = 0; i < tmpParticleNumber; i++) {
                tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
                tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
                tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            }

            // new ParticleIndexPairNumberCalculatorForTesting -------------------------------

            tmpStart = System.currentTimeMillis();
            ParticleIndexPairNumberCalculatorForTesting tmpParticleIndexPairNumberCalculator = new ParticleIndexPairNumberCalculatorForTesting(new Factory(Factory.RandomType.ACRNG_MT_64, 1, Factory.DpdType.CUTOFF_LENGTH_ONE, Factory.ElectrostaticsType.AD_HOC, Factory.BondType.HARMONIC, Factory.IntegrationType.GWMVV, new Double[] {0.65}), tmpMemoryLogger, tmpBoxSize, tmpPeriodicBoundaries, tmpCutOffLength, null, null, new ParallelizationInfo(2, 2, tmpParallelTaskNumber), tmpRandomNumberSeed);
            tmpEnd = System.currentTimeMillis();
            System.out.println(tmpOffset + "new tmpParticleIndexPairNumberCalculator,             t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));

            // Parallel cell-based ---------------------------------------------

            TestObjects tmpTestObjects = new TestObjects(
                new ConcurrentLinkedQueue<>(),
                new ConcurrentLinkedQueue<>(),
                new AtomicInteger(0),
                null
            );
            Parameters tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
            tmpStart = System.currentTimeMillis();
            tmpParticleIndexPairNumberCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
                tmpParticles.getR_x(), 
                tmpParticles.getR_y(), 
                tmpParticles.getR_z(), 
                tmpParameters
            );
            tmpEnd = System.currentTimeMillis();
            System.out.println(tmpOffset + "calculateCellBasedParticlePairInteractions(parallel), t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
            tmpNumber1 = tmpParameters.getTestObjects().getParticleIndexPairCounter().intValue();
            System.out.println(tmpOffset + "Number of particle index pairs                              = " + String.valueOf(tmpNumber1));
            System.out.println(tmpOffset + "Time in nanoseconds per particle index pair,         t [ns] = " + String.valueOf((int) ((double) (tmpEnd - tmpStart)/(double) tmpNumber1 * 1000000.0)));

            assertTrue("Parallel cell-based", tmpNumber1 > 0);

            System.out.println("");
        }        
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Test method
     * 
     * @param aParallelTaskNumber Parallel task number
     */
    private void test_IndexPairGeneration(int aParallelTaskNumber) {
        long tmpStart;
        long tmpEnd;
        String[] tmpParticleIndexPairStrings1;
        String[] tmpParticleIndexPairStrings2;

        System.out.println("----------------------------------------");
        System.out.println("test_IndexPairGeneration");
        System.out.println("aParallelTaskNumber = " + String.valueOf(aParallelTaskNumber));
        System.out.println("----------------------------------------");
        
        int tmpParticleNumber = 24000;
        AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
        double tmpBoxLength = 20.0;
        double tmpCutOffLength = 1.0;
        BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
        PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

        MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
        tmpMemoryLogger.start();
        ParticleArrays tmpParticles = new ParticleArrays(
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            null, null, null, null, null, null, null, null, null, null, null, null,
            null, null, null, null);
        ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
        for (int i = 0; i < tmpParticleNumber; i++) {
            tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
        }
        
        // new ParticleIndexPairCalculator -----------------------------------------
        
        tmpStart = System.currentTimeMillis();
        ParticleIndexPairCalculatorForTesting tmpParticleIndexPairCalculator = 
            new ParticleIndexPairCalculatorForTesting(
                new Factory(
                    Factory.RandomType.ACRNG_MT_64, 1, 
                    Factory.DpdType.CUTOFF_LENGTH_ONE, 
                    Factory.ElectrostaticsType.AD_HOC, 
                    Factory.BondType.HARMONIC, 
                    Factory.IntegrationType.GWMVV, 
                    new Double[] {0.65}
                ), 
                tmpMemoryLogger, 
                tmpBoxSize, 
                tmpPeriodicBoundaries, 
                tmpCutOffLength, 
                null, 
                null, 
                new ParallelizationInfo(2, 2, aParallelTaskNumber), 
                tmpRandomNumberSeed
            );
        tmpEnd = System.currentTimeMillis();
        System.out.println("new ParticleIndexPairCalculator,            t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        
        // Cell-based ----------------------------------------------------------

        TestObjects tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        Parameters tmpParameters = 
            new Parameters(
                null, 
                tmpParticles, 
                null, 
                null, 
                null, 
                null, 
                null, 
                tmpTestObjects
            );
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCellBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings1 = this.convertIntegerIndexPairQueue(tmpParameters.getTestObjects().getParticleIndexPairQueue());

        // Loop-based ----------------------------------------------------------

        tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        tmpParameters = 
            new Parameters(
                null, 
                tmpParticles, 
                null, 
                null, 
                null, 
                null, 
                null, 
                tmpTestObjects
            );
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairCalculator.calculateLoopBasedParticlePairInteractions(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateLoopBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings2 = this.convertIntegerIndexPairQueue(tmpParameters.getTestObjects().getParticleIndexPairQueue());
        
        assertTrue("Cell-based = Loop-based", this.isEqual(tmpParticleIndexPairStrings1, tmpParticleIndexPairStrings2));
        
        System.out.println("");
    }

    /**
     * Test method
     * 
     * @param aParallelTaskNumber Parallel task number
     */
    private void test_ParticlePairDistanceParametersCache(int aParallelTaskNumber) {
        long tmpStart;
        long tmpEnd;
        String[] tmpParticleIndexPairStrings1;
        String[] tmpParticleIndexPairStrings2;
        String[] tmpParticleIndexPairStrings1Cache;

        System.out.println("----------------------------------------");
        System.out.println("test_ParticlePairDistanceParametersCache");
        System.out.println("aParallelTaskNumber = " + String.valueOf(aParallelTaskNumber));
        System.out.println("----------------------------------------");
        
        int tmpParticleNumber = 24000;
        AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
        double tmpBoxLength = 20.0;
        double tmpCutOffLength = 1.0;
        BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
        PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

        MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
        ParticleArrays tmpParticles = new ParticleArrays(
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            null, null, null, null,
            null, null, null, null, null, null, null, null,
            null, null, null, null);
        ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
        for (int i = 0; i < tmpParticleNumber; i++) {
            tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
        }

        // NOTE: Integration type SCMVV is irrelevant, but cache is used in SCMVV
        Factory tmpFactory = 
            new Factory(
                Factory.RandomType.ACRNG_MT_64, 1, 
                Factory.DpdType.CUTOFF_LENGTH_ONE, 
                Factory.ElectrostaticsType.AD_HOC, 
                Factory.BondType.HARMONIC, 
                Factory.IntegrationType.SCMVV, 
                new Object[] {1, true}
            );
                
        // new ParticleIndexPairCalculator -------------------------------------
        
        tmpStart = System.currentTimeMillis();
        ParticleIndexPairCalculatorForTesting tmpParticleIndexPairCalculator = 
            new ParticleIndexPairCalculatorForTesting(
                tmpFactory, 
                tmpMemoryLogger, 
                tmpBoxSize, 
                tmpPeriodicBoundaries, 
                tmpCutOffLength, 
                null, 
                null, 
                new ParallelizationInfo(2, 2, aParallelTaskNumber), 
                tmpRandomNumberSeed
            );
        tmpEnd = System.currentTimeMillis();
        System.out.println("new ParticleIndexPairCalculator,             t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        
        // Sequential cell-based -----------------------------------------------

        TestObjects tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        Parameters tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairCalculator.setParticlePairDistanceParametersCacheActivity(true);
        tmpParticleIndexPairCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCellBasedParticlePairInteractions   t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings1 = this.convertIntegerIndexPairQueue(tmpParameters.getTestObjects().getParticleIndexPairQueue());
        ParticlePairDistanceParameters[][] tmpParticlePairDistanceParametersCache = tmpParticleIndexPairCalculator.getParticlePairDistanceParametersCache();
        tmpParticleIndexPairStrings1Cache = this.convertParticlePairDistanceParametersCache(tmpParticlePairDistanceParametersCache);

        assertTrue("Cell-based = Cache", this.isEqual(tmpParticleIndexPairStrings1, tmpParticleIndexPairStrings1Cache));

        // Cache-based ---------------------------------------------------------

        tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        ParticleIndexPairCalculatorForTesting tmpParticleIndexPairCalculatorWithCache = 
            new ParticleIndexPairCalculatorForTesting(
                tmpFactory, 
                tmpMemoryLogger, 
                tmpBoxSize, 
                tmpPeriodicBoundaries, 
                tmpCutOffLength, 
                null, 
                null, 
                new ParallelizationInfo(2, 2, aParallelTaskNumber), 
                tmpRandomNumberSeed
            );
        tmpParticleIndexPairCalculatorWithCache.setParticlePairDistanceParametersCache(tmpParticleIndexPairCalculator);
        
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairCalculatorWithCache.calculateCacheBasedParticlePairInteractions(tmpParameters);
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCacheBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings2 = this.convertIntegerIndexPairQueue(tmpParameters.getTestObjects().getParticleIndexPairQueue());
        
        assertTrue("Cell-based = Cache-based", this.isEqual(tmpParticleIndexPairStrings1, tmpParticleIndexPairStrings2));
        
        System.out.println("");
    }

    /**
     * Test method
     * 
     * @param aParallelTaskNumber Parallel task number
     */
    private void test_IndexPairRandomGeneration(int aParallelTaskNumber) {
        long tmpStart;
        long tmpEnd;
        String[] tmpParticleIndexPairStrings1;
        String[] tmpParticleIndexPairStrings2;

        System.out.println("----------------------------------------");
        System.out.println("test_IndexPairRandomGeneration");
        System.out.println("aParallelTaskNumber = " + String.valueOf(aParallelTaskNumber));
        System.out.println("----------------------------------------");
        
        int tmpParticleNumber = 24000;
        double tmpBoxLength = 20.0;
        double tmpCutOffLength = 1.0;
        BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
        PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

        MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
        ParticleArrays tmpParticles = new ParticleArrays(
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            null, null, null, null,
            null, null, null, null, null, null, null, null,
            null, null, null, null);
        ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
        for (int i = 0; i < tmpParticleNumber; i++) {
            tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
        }
        
        // new ParticleIndexPairCalculatorForTesting -----------------------------------------
        
        AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
        tmpStart = System.currentTimeMillis();
        ParticleIndexPairRandomCalculatorForTesting tmpParticleIndexPairRandomCalculator = new ParticleIndexPairRandomCalculatorForTesting(new Factory(Factory.RandomType.ACRNG_MT_64, 1, Factory.DpdType.CUTOFF_LENGTH_ONE, Factory.ElectrostaticsType.AD_HOC, Factory.BondType.HARMONIC, Factory.IntegrationType.GWMVV, new Double[] {0.65}), tmpMemoryLogger, tmpBoxSize, tmpPeriodicBoundaries, tmpCutOffLength, null, null, new ParallelizationInfo(2, 2, aParallelTaskNumber), tmpRandomNumberSeed);
        tmpEnd = System.currentTimeMillis();
        System.out.println("new tmpParticleIndexPairRandomCalculator,   t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        
        // Cell-based ----------------------------------------------------------

        TestObjects tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        Parameters tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairRandomCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCellBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings1 = tmpParameters.getTestObjects().getParticleIndexPairRandomQueue().toArray(new String[0]);
        Arrays.sort(tmpParticleIndexPairStrings1);

        // new ParticleIndexPairCalculatorForTesting -----------------------------------------
        
        tmpRandomNumberSeed = new AtomicInteger(1);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairRandomCalculator = new ParticleIndexPairRandomCalculatorForTesting(new Factory(Factory.RandomType.ACRNG_MT_64, 1, Factory.DpdType.CUTOFF_LENGTH_ONE, Factory.ElectrostaticsType.AD_HOC, Factory.BondType.HARMONIC, Factory.IntegrationType.GWMVV, new Double[] {0.65}), tmpMemoryLogger, tmpBoxSize, tmpPeriodicBoundaries, tmpCutOffLength, null, null, new ParallelizationInfo(2, 2, aParallelTaskNumber), tmpRandomNumberSeed);
        tmpEnd = System.currentTimeMillis();
        System.out.println("new tmpParticleIndexPairRandomCalculator,   t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        
        // Cell-based ----------------------------------------------------------

        tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairRandomCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCellBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpParticleIndexPairStrings2 = tmpParameters.getTestObjects().getParticleIndexPairRandomQueue().toArray(new String[0]);
        Arrays.sort(tmpParticleIndexPairStrings2);
        
        assertTrue("First = Second", this.isEqual(tmpParticleIndexPairStrings1, tmpParticleIndexPairStrings2));
        
        System.out.println("");
    }

    /**
     * Test method
     * 
     * @param aParallelTaskNumber Parallel task number
     */
    private void test_IndexPairCount(int aParallelTaskNumber) {
        long tmpStart;
        long tmpEnd;
        int tmpNumber1;
        int tmpNumber2;

        System.out.println("----------------------------------------");
        System.out.println("test_IndexPairCount");
        System.out.println("aParallelTaskNumber = " + String.valueOf(aParallelTaskNumber));
        System.out.println("----------------------------------------");
        
        int tmpParticleNumber = 24000;
        AtomicInteger tmpRandomNumberSeed = new AtomicInteger(1);
        double tmpBoxLength = 20.0;
        double tmpCutOffLength = 1.0;
        BoxSize tmpBoxSize = new BoxSize(tmpBoxLength, tmpBoxLength, tmpBoxLength);
        PeriodicBoundaries tmpPeriodicBoundaries = new PeriodicBoundaries(true, true, true);

        MemoryLogger tmpMemoryLogger = new MemoryLogger(new int[] {ILogger.EXCEPTION});
        ParticleArrays tmpParticles = new ParticleArrays(
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            new double[tmpParticleNumber],
            null, null, null, null,
            null, null, null, null, null, null, null, null,
            null, null, null, null);
        ApacheCommonsRandom tmpRandomDpd = new ApacheCommonsRandom(RandomSource.JDK, 1, 1);
        for (int i = 0; i < tmpParticleNumber; i++) {
            tmpParticles.getR_x()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_y()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
            tmpParticles.getR_z()[i] = tmpBoxLength * tmpRandomDpd.nextDouble();
        }
        
        // new ParticleIndexPairNumberCalculatorForTesting ---------------------
        
        tmpStart = System.currentTimeMillis();
        ParticleIndexPairNumberCalculatorForTesting tmpParticleIndexPairNumberCalculator = new ParticleIndexPairNumberCalculatorForTesting(new Factory(Factory.RandomType.ACRNG_MT_64, 1, Factory.DpdType.CUTOFF_LENGTH_ONE, Factory.ElectrostaticsType.AD_HOC, Factory.BondType.HARMONIC, Factory.IntegrationType.GWMVV, new Double[] {0.65}), tmpMemoryLogger, tmpBoxSize, tmpPeriodicBoundaries, tmpCutOffLength, null, null, new ParallelizationInfo(2, 2, aParallelTaskNumber), tmpRandomNumberSeed);
        tmpEnd = System.currentTimeMillis();
        System.out.println("new ParticleIndexPairNumberCalculator,      t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));

        // Sequential cell-based -----------------------------------------------

        TestObjects tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        Parameters tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairNumberCalculator.calculateCellBasedParticlePairInteractionsWithParticleCellAssignments(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateCellBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpNumber1 = tmpParameters.getTestObjects().getParticleIndexPairCounter().intValue();

        // Sequential loop -----------------------------------------------------

        tmpTestObjects = new TestObjects(
            new ConcurrentLinkedQueue<>(),
            new ConcurrentLinkedQueue<>(),
            new AtomicInteger(0),
            null
        );
        tmpParameters = new Parameters(null, tmpParticles, null, null, null, null, null, tmpTestObjects);
        tmpStart = System.currentTimeMillis();
        tmpParticleIndexPairNumberCalculator.calculateLoopBasedParticlePairInteractions(
            tmpParticles.getR_x(), 
            tmpParticles.getR_y(), 
            tmpParticles.getR_z(), 
            tmpParameters
        );
        tmpEnd = System.currentTimeMillis();
        System.out.println("calculateLoopBasedParticlePairInteractions, t [s] = " + String.valueOf((tmpEnd - tmpStart)/1000.0));
        tmpNumber2 = tmpParameters.getTestObjects().getParticleIndexPairCounter().intValue();
        
        assertTrue("loop = cell-based", tmpNumber2 == tmpNumber1);
        
        System.out.println("");
    }
    
    /**
     * Converts integer index pairs to corresponding strings where the 
     * smaller index is always the first in the string. The returned array is 
     * sorted ascending.
     * 
     * @param aParticleIndexPairQueue Queue with integer index pairs
     * @return Sorted integer index pair string array
     */
    private String[] convertIntegerIndexPairQueue(ConcurrentLinkedQueue<int[]> aParticleIndexPairQueue) {
        int[][] anIntegerIndexPairArray = aParticleIndexPairQueue.toArray(new int[0][]);
        String[] tmpIndexPairStringArray = new String[anIntegerIndexPairArray.length];
        for (int i = 0; i < anIntegerIndexPairArray.length; i++) {
            int[] tmpIndexPair = anIntegerIndexPairArray[i];
            if (tmpIndexPair[0] < tmpIndexPair[1]) {
                tmpIndexPairStringArray[i] = String.format("(%d, %d)", tmpIndexPair[0], tmpIndexPair[1]);
            } else {
                tmpIndexPairStringArray[i] = String.format("(%d, %d)", tmpIndexPair[1], tmpIndexPair[0]);
            }
        }
        Arrays.sort(tmpIndexPairStringArray);
        return tmpIndexPairStringArray;
    }

    /**
     * Converts integer index pairs in cache to corresponding strings where the 
     * smaller index is always the first in the string. The returned array is 
     * sorted ascending.
     * 
     * @param aParticlePairDistanceParametersCache ParticlePairDistanceParametersCache instance
     * @return Sorted integer index pair string array
     */
    private String[] convertParticlePairDistanceParametersCache(ParticlePairDistanceParameters[][] aParticlePairDistanceParametersCache) {
        LinkedList<String> tmpIndexPairStringList = new LinkedList<>();
        for (ParticlePairDistanceParameters[] tmpParticlePairDistanceParametersCachePerCellChunk : aParticlePairDistanceParametersCache) {
            for (ParticlePairDistanceParameters tmpParticlePairDistanceParameters : tmpParticlePairDistanceParametersCachePerCellChunk) {
                int[] tmpParticleIndex1 = tmpParticlePairDistanceParameters.getParticleIndices_i();
                int[] tmpParticleIndex2 = tmpParticlePairDistanceParameters.getParticleIndices_j();
                for (int k = 0; k < tmpParticlePairDistanceParameters.getSize(); k++) {
                    if (tmpParticleIndex1[k] < tmpParticleIndex2[k]) {
                        tmpIndexPairStringList.add(String.format("(%d, %d)", tmpParticleIndex1[k],  tmpParticleIndex2[k]));
                    } else {
                        tmpIndexPairStringList.add(String.format("(%d, %d)",  tmpParticleIndex2[k], tmpParticleIndex1[k]));
                    }
                }
            }
        }
        String[] tmpIndexPairStringArray = tmpIndexPairStringList.toArray(new String[0]);
        Arrays.sort(tmpIndexPairStringArray);
        return tmpIndexPairStringArray;
    }
    
    /**
     * Compares two string arrays
     * @param aStringArray1 String array 1
     * @param aStringArray2 String array 1
     * @return True: Both arrays are equal, false: Otherwise
     */
    private boolean isEqual(String[] aStringArray1, String[] aStringArray2) {
        if (aStringArray1.length != aStringArray2.length) {
            return false;
        }
        for (int i = 0; i < aStringArray1.length; i++) {
            if (!aStringArray1[i].equals(aStringArray2[i])) {
                return false;
            }
        }
        return true;
    }
     // </editor-fold>
   
}
