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

import de.gnwi.jdpd.interfaces.ILogLevel;
import de.gnwi.jdpdsp.samples.FileInput;
import de.gnwi.jdpdsp.samples.logger.FileLogger;
import de.gnwi.jdpdsp.samples.FileOutput;
import de.gnwi.jdpdsp.DpdSimulationTaskSP;
import de.gnwi.jdpdsp.samples.ProgressMonitor;
import junit.framework.TestCase;
import java.io.File;
import de.gnwi.jdpd.parameters.ParallelizationInfo;
import de.gnwi.jdpdsp.utilities.Utils;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Test class for DpdSimulationTask
 * 
 * @author Achim Zielesny
 */
public class DpdSimulationTaskTest extends TestCase {
    
    // <editor-fold defaultstate="collapsed" desc="Public tests">
    /**
     * Test method
     */
    public void test_Water_GWMVV() {
        this.test_Simulation("P:\\MFsim\\Code\\JdpdSP\\test\\de\\gnwi\\jdpdsp\\test_Water_GWMVV");
    }

    /**
     * Test method
     */
    public void test_IonPairs() {
        this.test_Simulation("P:\\MFsim\\Code\\JdpdSP\\test\\de\\gnwi\\jdpdsp\\test_IonPairs");
    }
    
    /**
     * Test method
     */
    public void test_DMPC() {
        this.test_Simulation("P:\\MFsim\\Code\\JdpdSP\\test\\de\\gnwi\\jdpdsp\\test_DMPC");
    }
    
    /**
     * Test method
     */
    public void test_KalataB1() {
        this.test_Simulation("P:\\MFsim\\Code\\JdpdSP\\test\\de\\gnwi\\jdpdsp\\test_KalataB1");
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Test method
     * 
     * @param aDirectoryPath Directory path
     */
    private void test_Simulation(String aDirectoryPath) {
        try{
            ExecutorService tmpExecutorService = Executors.newSingleThreadExecutor();
            String tmpOutputDirectoryPath = aDirectoryPath + File.separatorChar + "output";
            String tmpLogFilePathname = tmpOutputDirectoryPath + File.separatorChar + "LogFile.txt";
            // Copy old log file to temporary directory
            String tmpOldLogFileName = "LogFile.txt.old";
            Utils.copySingleFileToTemporaryDirectory(tmpLogFilePathname, tmpOldLogFileName);
            Utils.clearDirectory(tmpOutputDirectoryPath);
            // File logger
            if ((new File(tmpLogFilePathname)).exists()) {
                (new File(tmpLogFilePathname)).delete();
            }
            int[] tmpLogLevels = ILogLevel.ALL_LOGLEVELS;
            // int[] tmpLogLevels = 
            //     new int[] {
            //         ILogger.EXCEPTION, 
            //         ILogger.A_IJ, 
            //         ILogger.OUTPUT_STEP,
            //         ILogger.QUANTITY,
            //         ILogLevel.PARTICLE
            //     };
            FileLogger tmpFileLogger = new FileLogger(tmpLogFilePathname, tmpLogLevels);
            // File input
            String tmpInputFilePathname = aDirectoryPath + File.separatorChar + "Input.txt";
            FileInput tmpFileInput = new FileInput(tmpInputFilePathname, tmpFileLogger);
            // File output
            FileOutput tmpFileOutput = 
                new FileOutput(
                    tmpOutputDirectoryPath,
                    tmpOutputDirectoryPath,
                    null,
                    null,
                    null,
                    null,
                    16,
                    -1
                );
            // Progress monitor
            ProgressMonitor tmpProgressMonitor = new ProgressMonitor();
            // DpdSimulationTask
            DpdSimulationTaskSP tmpDpdSimulationTask = new DpdSimulationTaskSP(null, tmpFileInput, tmpFileOutput, tmpProgressMonitor, tmpFileLogger, new ParallelizationInfo(100, 100, 8));
            Future<Boolean> tmpFuture = tmpExecutorService.submit(tmpDpdSimulationTask);
            while (!tmpProgressMonitor.hasFinished()) {
                Thread.sleep(100);
            }
            // Finish
            tmpExecutorService.shutdown();
            // Move old log file from temporary directory
            Utils.moveSingleFileFromTemporaryDirectory(tmpOldLogFileName, tmpOutputDirectoryPath);
            assertTrue("NoException", true);
        } catch (Exception anException) {
            assertTrue("Exception", false);
        }
    }
    // </editor-fold>

}
