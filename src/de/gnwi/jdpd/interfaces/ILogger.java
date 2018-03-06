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
package de.gnwi.jdpd.interfaces;

/**
 * Interface for simulation logger
 * 
 * @author Achim Zielesny
 */
public interface ILogger {

    // <editor-fold defaultstate="collapsed" desc="Enum LogLevel">
    /**
     * Log level
     */
    public enum LogLevel {
            
        /**
         * Log level EXCEPTION: Exceptions are logged
         */
        EXCEPTIONS,
        /**
         * Log level SIMULATION_INIT: Simulation initialisation is logged
         */
        SIMULATION_INIT,
        /**
         * Log level SIMULATION_PROGRESS: Simulation progress is logged
         */
        SIMULATION_PROGRESS,
        /**
         * Log level TIME_STEPS: Time steps are logged
         */
        TIME_STEPS,
        /**
         * Log level OUTPUT_TIME_STEPS: Output time steps are logged
         */
        OUTPUT_TIME_STEPS,
        /**
         * Log level METHOD_CALLS: Each method entry/exit is logged
         */
        METHOD_CALLS,
        /**
         * Log level VELOCITY_SCALE_FACTOR: Velocity scale factor is logged.
         */
        VELOCITY_SCALE_FACTOR,
        /**
         * Log level INTERMEDIATE_RESULTS: Intermediate results are logged.
         */
        INTERMEDIATE_RESULTS,
        /**
         * Log level PARALLELIZATION: Parallelisation related quantities are logged.
         */
        PARALLELIZATION,
        /**
         * Log level A_IJ: a(ij) are logged.
         */
        A_IJ

    }    
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Checks if current log levels comprise aLogLevel
     * NOTE: A fast method implementation is MANDATORY.
     * 
     * @param aLogLevel Log level to be tested
     * @return True: Current log level comprises aLogLevel, false: Otherwise.
     */
    boolean isLogLevel(LogLevel aLogLevel);

    /**
     * Starts logger
     */
    void start();

    // <editor-fold defaultstate="collapsed" desc="LogLevel EXCEPTIONS related methods">
    /**
     * Appends with log level EXCEPTIONS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMethod  Method in which exception occurred
     * @param aStacktrace Stack trace of excpetion
     */
    public void appendException(String aMethod, String aStacktrace);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel SIMULATION_INIT related methods">
    /**
     * Appends with log level SIMULATION_INIT
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationInit(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel SIMULATION_PROGRESS related methods">
    /**
     * Appends with log level SIMULATION_PROGRESS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationProgress(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel TIME_STEPS related methods">
    /**
     * Appends with log level TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendTimeStep(String aMessage);

    /**
     * Appends start of time step with log level TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepStart(int aTimeStep, long anId);

    /**
     * Appends end of time step with log level TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepEnd(int aTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel OUTPUT_TIME_STEPS related methods">
    /**
     * Appends with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendOutputTimeStep(String aMessage);

    /**
     * Appends start of output time step with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepStart(int anOutputTimeStep, long anId);

    /**
     * Appends end of output time step with log level OUTPUT_TIME_STEPS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepEnd(int anOutputTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel METHOD_CALLS related methods">
    /**
     * Appends with log level METHOD_CALLS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendMethodCall(String aMessage);

    /**
     * Appends start of method with log level METHOD_CALLS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallStart(String aMessage, long anId);

    /**
     * Appends end of method with log level METHOD_CALLS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallEnd(String aMessage, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel VELOCITY_SCALE_FACTOR related methods">
    /**
     * Appends with log level VELOCITY_SCALE_FACTOR
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendVelocityScaleFactor(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel INTERMEDIATE_RESULTS related methods">
    /**
     * Appends with log level INTERMEDIATE_RESULTS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendIntermediateResults(String aMessage);
    /**
     * Appends with log level INTERMEDIATE_RESULTS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anArrayName Array name
     * @param anArray Array
     */
    void appendIntermediateResults(String anArrayName, double[] anArray);    
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel PARALLELIZATION related methods">
    /**
     * Appends with log level PARALLELIZATION
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendParallelization(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel A_IJ related methods">
    /**
     * Appends with log level A_IJ
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendAij(String aMessage);
    // </editor-fold>

    /**
     * Returns ID
     * 
     * @return ID
     */
    long getId();
    
    /**
     * Finishes logger
     * 
     * @return true: Operation was successful, false: Operation failed
     */
    boolean finish();
    // </editor-fold>
    
}
