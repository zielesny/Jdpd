/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpdsp.samples.logger;

import de.gnwi.jdpd.interfaces.ILogLevel;
import java.util.concurrent.ConcurrentLinkedQueue;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.utilities.RegexPatterns;
import de.gnwi.jdpdsp.utilities.Utils;
import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.DoubleAdder;
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
    public final String FORMAT = "%s : %s";

    /**
     * Logger format string with ID
     */
    public final String ID_FORMAT = "%s %s : %s (ID = %s)";
    
    /**
     * Start prefix for simulation logger
     */
    public final String START_PREFIX = "<";
    
    /**
     * End prefix for simulation logger
     */
    public final String END_PREFIX = ">";
    
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

    /**
     * True: SCMVV information accumulation, false: Otherwise
     */
    private boolean isScmvvInformationAccumulation;
    
    /**
     * SCMVV dissipative force particle index pair counter
     */
    private AtomicInteger scmvvFdissParticleIndexPairCounter;

    /**
     * SCMVV dissipative force Rij_x adder
     */
    private DoubleAdder scmvvFdissRij_x_Adder;

    /**
     * SCMVV dissipative force Abs(Rij_x) adder
     */
    private DoubleAdder scmvvFdissAbsRij_x_Adder;

    /**
     * SCMVV dissipative force Vij_x adder
     */
    private DoubleAdder scmvvFdissVij_x_Adder;

    /**
     * SCMVV dissipative force Abs(Vij_x) adder
     */
    private DoubleAdder scmvvFdissAbsVij_x_Adder;

    /**
     * SCMVV dissipative force Rij_x*Vij_x adder
     */
    private DoubleAdder scmvvFdissRijVij_x_Adder;

    /**
     * SCMVV dissipative force Abs(Rij_x*Vij)_x adder
     */
    private DoubleAdder scmvvFdissAbsRijVij_x_Adder;

    /**
     * SCMVV dissipative force VdotR adder
     */
    private DoubleAdder scmvvFdissVdotR_Adder;

    /**
     * SCMVV dissipative force gamma factor adder
     */
    private DoubleAdder scmvvFdissGammaFactor_Adder;
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
    protected boolean[] logLevelArray;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aLogLevels Log levels
     * @throws IllegalArgumentException Thrown if argument is illegal
     */
    public MemoryLogger(int[] aLogLevels) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aLogLevels == null || aLogLevels.length == 0) {
            throw new IllegalArgumentException("MemoryLogger.Constructor: aLogLevels is null/length 0.");
        }
        // </editor-fold>
        this.logLevelArray = new boolean[ILogLevel.LOGLEVEL_COUNT];
        Arrays.fill(this.logLevelArray, false);
        for (int tmpLogLevel : aLogLevels) {
            this.logLevelArray[tmpLogLevel] = true;
        }
        this.indentingOffset = "";
        this.id = new AtomicLong();
        this.logQueue = new ConcurrentLinkedQueue<>();
        
        // Log level SCMVV related class variables
        this.isScmvvInformationAccumulation = false;
        this.scmvvFdissParticleIndexPairCounter = null;
        this.scmvvFdissRij_x_Adder = null;
        this.scmvvFdissAbsRij_x_Adder = null;
        this.scmvvFdissVij_x_Adder = null;
        this.scmvvFdissAbsVij_x_Adder = null;
        this.scmvvFdissRijVij_x_Adder = null;
        this.scmvvFdissAbsRijVij_x_Adder = null;
        this.scmvvFdissVdotR_Adder = null;
        this.scmvvFdissGammaFactor_Adder = null;
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
    public boolean isLogLevel(int aLogLevel) {
        return this.logLevelArray[aLogLevel];
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

    // <editor-fold defaultstate="collapsed" desc="- LogLevel EXCEPTION related methods">
    /**
     * Appends with log level EXCEPTION
     * NOTE: No checks are performed
     * 
     * @param aMethod  Method in which exception occurred
     * @param aStacktrace Stack trace of excpetion
     */
    @Override
    public void appendException(String aMethod, String aStacktrace) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.EXCEPTION)) {
            this.reentrantLock.lock();
            try {
                long tmpId = this.getId();
                StringBuilder tmpBuffer = new StringBuilder(aStacktrace.length() + aMethod.length() + 100);
                tmpBuffer.append(String.format(this.ID_FORMAT, this.START_PREFIX, "EXCEPTION", aMethod, String.valueOf(tmpId)));
                tmpBuffer.append("\n");
                tmpBuffer.append(aStacktrace);
                tmpBuffer.append(String.format(this.ID_FORMAT, this.END_PREFIX, "EXCEPTION", aMethod, String.valueOf(tmpId)));
                this.appendException(tmpBuffer.toString());
            } finally {
                this.reentrantLock.unlock();
            }
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel INIT related methods">
    /**
     * Appends with log level INIT
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendSimulationInit(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.INIT)) {
            this.append(String.format(this.FORMAT, "INIT       ", aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel PROGRESS related methods">
    /**
     * Appends with log level PROGRESS
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendSimulationProgress(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.PROGRESS)) {
            this.append(String.format(this.FORMAT, "PROGRESS   ", aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel TIME_STEP related methods">
    /**
     * Appends with log level TIME_STEP
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendTimeStep(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.TIME_STEP)) {
            this.append(String.format(this.FORMAT, "TIME_STEP  ", aMessage));
        }
    }

    /**
     * Appends start of time step with log level TIME_STEP
     * NOTE: No checks are performed
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    @Override
    public void appendTimeStepStart(int aTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.TIME_STEP)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, "TIME_STEP  ", "Time step = " + String.valueOf(aTimeStep), String.valueOf(anId)));
        }
    }

    /**
     * Appends end of time step with log level TIME_STEP
     * NOTE: No checks are performed
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    @Override
    public void appendTimeStepEnd(int aTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.TIME_STEP)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, "TIME_STEP  ", "Time step = " + String.valueOf(aTimeStep), String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel OUTPUT_STEP related methods">
    /**
     * Appends with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendOutputTimeStep(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.OUTPUT_STEP)) {
            this.append(String.format(this.FORMAT, "OUTPUT_STEP", aMessage));
        }
    }

    /**
     * Appends start of output time step with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    @Override
    public void appendOutputTimeStepStart(int anOutputTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.OUTPUT_STEP)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, "OUTPUT_STEP", "Output time step = " + String.valueOf(anOutputTimeStep), String.valueOf(anId)));
        }
    }

    /**
     * Appends end of output time step with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    @Override
    public void appendOutputTimeStepEnd(int anOutputTimeStep, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.OUTPUT_STEP)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, "OUTPUT_STEP", "Output time step = " + String.valueOf(anOutputTimeStep), String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel METHOD_CALL related methods">
    /**
     * Appends with log level METHOD_CALL
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     */
    @Override
    public void appendMethodCall(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.METHOD_CALL)) {
            this.append(String.format(this.FORMAT, "METHOD_CALL", aMessage));
        }
    }

    /**
     * Appends start of method with log level METHOD_CALL
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    @Override
    public void appendMethodCallStart(String aMessage, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.METHOD_CALL)) {
            this.append(String.format(this.ID_FORMAT, this.START_PREFIX, "METHOD_CALL", aMessage, String.valueOf(anId)));
        }
    }

    /**
     * Appends end of method with log level METHOD_CALL
     * NOTE: No checks are performed
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    @Override
    public void appendMethodCallEnd(String aMessage, long anId) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.METHOD_CALL)) {
            this.append(String.format(this.ID_FORMAT, this.END_PREFIX, "METHOD_CALL", aMessage, String.valueOf(anId)));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel V_SCALE related methods">
    /**
     * Appends with log level V_SCALE
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendVelocityScaleFactor(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.V_SCALE)) {
            this.append(String.format(this.FORMAT, "V_SCALE    ", aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel QUANTITY related methods">
    /**
     * Appends with log level QUANTITY
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendIntermediateResults(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.QUANTITY)) {
            this.append(String.format(this.FORMAT, "QUANTITY   ", aMessage));
        }
    }
    /**
     * Appends with log level QUANTITY
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anArrayName Array name
     * @param anArray Double array
     */
    @Override
    public void appendIntermediateResults(String anArrayName, float[] anArray) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.QUANTITY)) {
            String tmpFormat = anArrayName + "[%s] = %s";
            for (int i = 0; i < anArray.length; i++) {
                String tmpMessage = String.format(tmpFormat, String.valueOf(i), String.valueOf(anArray[i]));
                this.append(String.format(this.FORMAT, "QUANTITY   ", tmpMessage));
            }
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel PARALLEL related methods">
    /**
     * Appends with log level PARALLEL
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendParallelization(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.PARALLEL)) {
            this.append(String.format(this.FORMAT, "PARALLEL   ", aMessage));
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
        if (this.isLogLevel(ILogLevel.A_IJ)) {
            this.append(String.format(this.FORMAT, "A_IJ       ", aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel PARTICLE related methods">
    /**
     * Appends with log level PARTICLE
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendParticleForceMagnitude(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.PARTICLE)) {
            this.append(String.format(this.FORMAT, "PARTICLE   ", aMessage));
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="- LogLevel SCMVV related methods">
    /**
     * Appends with log level SCMVV
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    @Override
    public void appendScmvv(String aMessage) {
        // No checks are performed!
        if (this.isLogLevel(ILogLevel.SCMVV)) {
            this.append(String.format(this.FORMAT, "SCMVV", aMessage));
        }
    }

    /**
     * True: SCMVV information accumulation, false: Otherwise
      * 
     * @param aValue Value
     */
    @Override
    public void setScmvvInformationAccumulation(boolean aValue) {
        this.isScmvvInformationAccumulation = aValue;
    }
    
    /**
     * True: SCMVV information accumulation, false: Otherwise
      * 
     * @return True: SCMVV information accumulation, false: Otherwise
     */
    @Override
    public boolean isScmvvInformationAccumulation() {
        return this.isScmvvInformationAccumulation;
    }
    
    /**
     * SCMVV dissipative force particle index pair counter
     * 
     * @return SCMVV dissipative force particle index pair counter
     */
    @Override
    public AtomicInteger getScmvvFdissParticleIndexPairCounter() {
        return this.scmvvFdissParticleIndexPairCounter;
    }
    
    /**
     * SCMVV dissipative force particle index pair counter
     * 
     * @param aParticleIndexPairCounter Particle index pair counter
     */
    @Override
    public void setScmvvFdissParticleIndexPairCounter(AtomicInteger aParticleIndexPairCounter) {
        this.scmvvFdissParticleIndexPairCounter = aParticleIndexPairCounter;
    }

    /**
     * SCMVV dissipative force Rij_x adder
     * 
     * @return SCMVV dissipative force Rij_x adder
     */
    @Override
    public DoubleAdder getScmvvFdissRij_x_Adder() {
        return this.scmvvFdissRij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Rij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissRij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissRij_x_Adder = aDoubleAdder;
    }

    /**
     * SCMVV dissipative force Abs(Rij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Rij_x) adder
     */
    @Override
    public DoubleAdder getScmvvFdissAbsRij_x_Adder() {
        return this.scmvvFdissAbsRij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Abs(Rij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissAbsRij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissAbsRij_x_Adder = aDoubleAdder;
    }

    /**
     * SCMVV dissipative force Vij_x adder
     * 
     * @return SCMVV dissipative force Vij_x adder
     */
    @Override
    public DoubleAdder getScmvvFdissVij_x_Adder() {
        return this.scmvvFdissVij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Vij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissVij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissVij_x_Adder = aDoubleAdder;
    }
    
    /**
     * SCMVV dissipative force Abs(Vij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Vij_x) adder
     */
    @Override
    public DoubleAdder getScmvvFdissAbsVij_x_Adder() {
        return this.scmvvFdissAbsVij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Abs(Vij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissAbsVij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissAbsVij_x_Adder = aDoubleAdder;
    }

    /**
     * SCMVV dissipative force Rij_x*Vij_x adder
     * 
     * @return SCMVV dissipative force Rij_x*Vij_x adder
     */
    @Override
    public DoubleAdder getScmvvFdissRijVij_x_Adder() {
        return this.scmvvFdissRijVij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Rij_x*Vij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissRijVij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissRijVij_x_Adder = aDoubleAdder;
    }

    /**
     * SCMVV dissipative force Abs(Rij_x*Vij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Rij_x*Vij)_x adder
     */
    @Override
    public DoubleAdder getScmvvFdissAbsRijVij_x_Adder() {
        return this.scmvvFdissAbsRijVij_x_Adder;
    }
    
    /**
     * SCMVV dissipative force Abs(Rij_x*Vij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissAbsRijVij_x_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissAbsRijVij_x_Adder = aDoubleAdder;
    }
    
    /**
     * SCMVV dissipative force VdotR adder
     * 
     * @return SCMVV dissipative force VdotR adder
     */
    @Override
    public DoubleAdder getScmvvFdissVdotR_Adder() {
        return this.scmvvFdissVdotR_Adder;
    }
    
    /**
     * SCMVV dissipative force VdotR adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissVdotR_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissVdotR_Adder = aDoubleAdder;
    }

    /**
     * SCMVV dissipative force GammaFactor adder
     * 
     * @return SCMVV dissipative force GammaFactor adder
     */
    @Override
    public DoubleAdder getScmvvFdissGammaFactor_Adder() {
        return this.scmvvFdissGammaFactor_Adder;
    }
    
    /**
     * SCMVV dissipative force GammaFactor adder
     * 
     * @param aDoubleAdder Double adder
     */
    @Override
    public void setScmvvFdissGammaFactor_Adder(DoubleAdder aDoubleAdder) {
        this.scmvvFdissGammaFactor_Adder = aDoubleAdder;
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
     * (No checks are performed)
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
     * (No checks are performed)
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
