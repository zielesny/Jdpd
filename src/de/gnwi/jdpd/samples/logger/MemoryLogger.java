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
package de.gnwi.jdpd.samples.logger;

import java.util.HashMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import de.gnwi.jdpd.interfaces.ILogger;
import de.gnwi.jdpd.utilities.RegexPatterns;
import de.gnwi.jdpd.utilities.Utils;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

/**
 * File logger: Write log information to file when finishing.
 * 
 * @author Achim Zielesny
 */
public class MemoryLogger implements ILogger {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    // <editor-fold defaultstate="collapsed" desc="- Atomics">
    /**
     * ID
     */
    private final AtomicLong id;
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Strings">
    /**
     * Logger format string
     */
    public final String FORMAT = "(%s): %s";

    /**
     * Logger format string with ID
     */
    public final String ID_FORMAT = "%s(%s): %s (ID = %s)";
    
    /**
     * Start prefix for simulation logger
     */
    public final String START_PREFIX = "<";
    
    /**
     * End prefix for simulation logger
     */
    public final String END_PREFIX = ">";
    
    /**
     * Abbreviation for exception
     */
    public final String ABBREVIATION_EXCEPTION = "EXCEPTION";
    
    /**
     * Abbreviation for simulation initialisation
     */
    public final String ABBREVIATION_SIMULATION_INIT = "SIM_INIT";
    
    /**
     * Abbreviation for time step
     */
    public final String ABBREVIATION_TIME_STEPS = "TIME_STEP";
    
    /**
     * Abbreviation for output time step
     */
    public final String ABBREVIATION_OUTPUT_TIME_STEPS = "OUTPUT_STEP";
    
    /**
     * Abbreviation for simulation progress
     */
    public final String ABBREVIATION_SIMULATION_PROGRESS = "SIM_PROGRESS";
    
    /**
     * Abbreviation for method call
     */
    public final String ABBREVIATION_METHOD_CALL = "METHOD";
    
    /**
     * Abbreviation for velocity scale factor
     */
    public final String ABBREVIATION_VELOCITY_SCALE_FACTOR = "SCALE_FACTOR";
    
    /**
     * Abbreviation for intermediate results
     */
    public final String ABBREVIATION_INTERMEDIATE_RESULTS = "RESULT";
    
    /**
     * Abbreviation for parallelisation
     */
    public final String ABBREVIATION_PARALLELIZATION = "PARALLEL";
    
    /**
     * Abbreviation for a(ij)
     */
    public final String ABBREVIATION_A_IJ = "A_IJ";
    
    /**
     * Logger indent string
     */
    public final String LOGGER_INDENT = "|  ";
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- Miscellaneous">
    private final Lock reentrantLock = new ReentrantLock();
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * String offset for indenting
     */
    private String indentingOffset;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected class variables">
    /**
     * Queue for logging
     */
    protected ConcurrentLinkedQueue<String> logQueue;
    
    /**
     * Hash map for log levels
     */
    protected HashMap<ILogger.LogLevel, ILogger.LogLevel> logLevelMap;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aLogLevels Log levels
     * @throws IllegalArgumentException Thrown if argument is illegal
     */
    public MemoryLogger(ILogger.LogLevel[] aLogLevels) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aLogLevels == null || aLogLevels.length == 0) {
            throw new IllegalArgumentException("MemoryLogger.Constructor: aLogLevels is null/length 0.");
        }
        // </editor-fold>
        this.logLevelMap = new HashMap<>(aLogLevels.length);
        for (ILogger.LogLevel tmpLogLevel : aLogLevels) {
            this.logLevelMap.put(tmpLogLevel, tmpLogLevel);
        }
        this.indentingOffset = "";
        this.id = new AtomicLong();
        this.logQueue = new ConcurrentLinkedQueue<>();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public overridden methods">
    /**
     * Checks if current log level comprises aLogLevel
     * NOTE: This method should NOT be synchronised for better performance.
     * 
     * @param aLogLevel Log level to be tested
     * @return True: Current log level comprises aLogLevel, false: Otherwise.
     */
    @Override
    public boolean isLogLevel(ILogger.LogLevel aLogLevel) {
        return this.logLevelMap.containsKey(aLogLevel);
    }

    /**
     * Starts logger
     */
    @Override
    public void start() {
        this.logQueue.clear();
        this.logQueue.add("---------------------------------------");
        this.logQueue.add("DPD Simulation:");
        this.logQueue.add("Logger started at " + Utils.getTimestamp());
        this.logQueue.add("---------------------------------------");
    }

    // <editor-fold defaultstate="collapsed" desc="- LogLevel EXCEPTIONS related methods">
    /**
     * Appends with log level EXCEPTIONS
     * NOTE: No checks are performed
     * 
     * @param aMethod  Method in which exception occurred
     * @param aStacktrace Stack trace of excpetion
     */
    @Override
    public void appendException(String aMethod, String aStacktrace) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.EXCEPTIONS)) {
            this.reentrantLock.lock();
            try {
                long tmpId = this.getId();
                StringBuilder tmpBuffer = new StringBuilder(aStacktrace.length() + aMethod.length() + 100);
                tmpBuffer.append(String.format(this.ID_FORMAT, this.START_PREFIX, this.ABBREVIATION_EXCEPTION, aMethod, String.valueOf(tmpId)));
                tmpBuffer.append("\n");
                tmpBuffer.append(aStacktrace);
                tmpBuffer.append(String.format(this.ID_FORMAT, this.END_PREFIX, this.ABBREVIATION_EXCEPTION, aMethod, String.valueOf(tmpId)));
                this.appendException(tmpBuffer.toString());
            } finally {
                this.reentrantLock.unlock();
            }
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel SIMULATION_INIT related methods">
    /**
     * Appends with log level SIMULATION_INIT
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendSimulationInit(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.SIMULATION_INIT)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_SIMULATION_INIT, aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel SIMULATION_PROGRESS related methods">
    /**
     * Appends with log level SIMULATION_PROGRESS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendSimulationProgress(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.SIMULATION_PROGRESS)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_SIMULATION_PROGRESS, aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel TIME_STEPS related methods">
    /**
     * Appends with log level TIME_STEPS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendTimeStep(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.TIME_STEPS)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_TIME_STEPS, aMessage));
        }
    }

    /**
     * Appends start of time step with log level TIME_STEPS
     * NOTE: No checks are performed
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    @Override
    public void appendTimeStepStart(int aTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.TIME_STEPS)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, this.ABBREVIATION_TIME_STEPS, "Time step = " + String.valueOf(aTimeStep), String.valueOf(anId)));
        }
    }

    /**
     * Appends end of time step with log level TIME_STEPS
     * NOTE: No checks are performed
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    @Override
    public void appendTimeStepEnd(int aTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.TIME_STEPS)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, this.ABBREVIATION_TIME_STEPS, "Time step = " + String.valueOf(aTimeStep), String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel OUTPUT_TIME_STEPS related methods">
    /**
     * Appends with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendOutputTimeStep(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.OUTPUT_TIME_STEPS)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_OUTPUT_TIME_STEPS, aMessage));
        }
    }

    /**
     * Appends start of output time step with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    @Override
    public void appendOutputTimeStepStart(int anOutputTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.OUTPUT_TIME_STEPS)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, this.ABBREVIATION_OUTPUT_TIME_STEPS, "Output time step = " + String.valueOf(anOutputTimeStep), String.valueOf(anId)));
        }
    }

    /**
     * Appends end of output time step with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    @Override
    public void appendOutputTimeStepEnd(int anOutputTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.OUTPUT_TIME_STEPS)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, this.ABBREVIATION_OUTPUT_TIME_STEPS, "Output time step = " + String.valueOf(anOutputTimeStep), String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel METHOD_CALLS related methods">
    /**
     * Appends with log level METHOD_CALLS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendMethodCall(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.METHOD_CALLS)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_METHOD_CALL, aMessage));
        }
    }

    /**
     * Appends start of method with log level METHOD_CALLS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    @Override
    public void appendMethodCallStart(String aMessage, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.METHOD_CALLS)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, this.ABBREVIATION_METHOD_CALL, aMessage, String.valueOf(anId)));
        }
    }

    /**
     * Appends end of method with log level METHOD_CALLS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    @Override
    public void appendMethodCallEnd(String aMessage, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.METHOD_CALLS)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, this.ABBREVIATION_METHOD_CALL, aMessage, String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel VELOCITY_SCALE_FACTOR related methods">
    /**
     * Appends with log level VELOCITY_SCALE_FACTOR
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendVelocityScaleFactor(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.VELOCITY_SCALE_FACTOR)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_VELOCITY_SCALE_FACTOR, aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel INTERMEDIATE_RESULTS related methods">
    /**
     * Appends with log level INTERMEDIATE_RESULTS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendIntermediateResults(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.INTERMEDIATE_RESULTS)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_INTERMEDIATE_RESULTS, aMessage));
        }
    }
    /**
     * Appends with log level INTERMEDIATE_RESULTS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anArrayName Array name
     * @param anArray Double array
     */
    @Override
    public void appendIntermediateResults(String anArrayName, double[] anArray) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.INTERMEDIATE_RESULTS)) {
            String tmpFormat = anArrayName + "[%s] = %s";
            for (int i = 0; i < anArray.length; i++) {
                String tmpMessage = String.format(tmpFormat, String.valueOf(i), String.valueOf(anArray[i]));
                this.append(String.format(this.FORMAT, this.ABBREVIATION_INTERMEDIATE_RESULTS, tmpMessage));
            }
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel PARALLELIZATION related methods">
    /**
     * Appends with log level PARALLELIZATION
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendParallelization(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.PARALLELIZATION)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_PARALLELIZATION, aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel A_IJ related methods">
    /**
     * Appends with log level A_IJ
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendAij(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogger.LogLevel.A_IJ)) {
            this.append(String.format(this.FORMAT, this.ABBREVIATION_A_IJ, aMessage));
        }
    }
    // </editor-fold>

    /**
     * Returns ID
     * 
     * @return ID
     */
    @Override
    public long getId() {
        return this.id.incrementAndGet();
    }
    
    /**
     * Finishes logger
     * 
     * @return true: Operation was successful, false: Operation failed
     */
    @Override
    public boolean finish() {
        this.logQueue.add("----------------------------------------");
        this.logQueue.add("Logger finished at " + Utils.getTimestamp());
        this.logQueue.add("----------------------------------------");
        return true;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Log queue
     * 
     * @return Log queue
     */
    public ConcurrentLinkedQueue<String> getLogQueue() {
        return this.logQueue;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private synchronized methods">
    /**
     * Appends a message. If message contains several lines the message is split
     * in corresponding lines with each line appended separately.
     * NOTE: Method MUST be synchronised in order to get atomic 
     *       append-operations.
     * NOTE: No checks are performed.
     * 
     * @param aMessage Message
     */
    private synchronized void append(String aMessage) {
        // No checks are performed
        if (aMessage.startsWith(this.END_PREFIX)) {
            if (!this.indentingOffset.isEmpty()) {
                this.indentingOffset = this.indentingOffset.substring(0, this.indentingOffset.length() - this.LOGGER_INDENT.length());
            }
        }
        this.logQueue.add(this.indentingOffset + aMessage);
        if (aMessage.startsWith(this.START_PREFIX)) {
            this.indentingOffset += this.LOGGER_INDENT;
        }
    }

    /**
     * Appends an exception.
     * NOTE: Method MUST be synchronised in order to get atomic 
     *       append-operations.
     * NOTE: No checks are performed.
     * 
     * @param aMessage Message
     */
    private synchronized void appendException(String aMessage) {
        // No checks are performed
        this.logQueue.add("---------------------------------------");
        String[] tmpMessageLines = RegexPatterns.NEW_LINE_PATTERN.split(aMessage);
        for (String tmpMessageLine : tmpMessageLines) {
            if (tmpMessageLine.startsWith(this.START_PREFIX) || tmpMessageLine.startsWith(this.END_PREFIX)) {
                this.logQueue.add(tmpMessageLine);
            } else {
                this.logQueue.add(this.LOGGER_INDENT + tmpMessageLine);
            }
        }
        this.logQueue.add("---------------------------------------");
        if (!this.indentingOffset.isEmpty()) {
            this.indentingOffset = this.indentingOffset.substring(0, this.indentingOffset.length() - this.LOGGER_INDENT.length());
        }
    }
    // </editor-fold>
    
}
