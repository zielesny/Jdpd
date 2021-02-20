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
package de.gnwi.jdpd.interfaces;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Interface for simulation logger
 * 
 * @author Achim Zielesny
 */
public interface ILogger {

    // <editor-fold defaultstate="collapsed" desc="Static final LogLevels">
    // IMPORTANT: Log levels MUST be consecutively numbered started with 0!
    /**
     * Log level EXCEPTION: Exceptions are logged
     */
    static final int EXCEPTION = 0;
    /**
     * Log level INIT: Simulation initialisation is logged
     */
    static final int INIT = 1;
    /**
     * Log level PROGRESS: Simulation progress is logged
     */
    static final int PROGRESS = 2;
    /**
     * Log level TIME_STEP: Time steps are logged
     */
    static final int TIME_STEP = 3;
    /**
     * Log level OUTPUT_STEP: Output time steps are logged
     */
    static final int OUTPUT_STEP = 4;
    /**
     * Log level METHOD_CALL: Each method entry/exit is logged
     */
    static final int METHOD_CALL = 5;
    /**
     * Log level V_SCALE: Velocity scale factor is logged.
     */
    static final int V_SCALE = 6;
    /**
     * Log level QUANTITY: Intermediate values of quantities are logged.
     */
    static final int QUANTITY = 7;
    /**
     * Log level PARALLEL: Parallelization related quantities are logged.
     */
    static final int PARALLEL = 8;
    /**
     * Log level A_IJ: a(ij) are logged.
     */
    static final int A_IJ = 9;
    /**
     * Log level PARTICLE: Particle properties (forces, velocities etc.) are analysed (min, mean, max)
     */
    static final int PARTICLE = 10;
    /**
     * Log level SCMVV: Specific logging for integration type SCMVV
     */
    static final int SCMVV = 11;
    /**
     * Array with all log levels
     * NOTE: MUST be correct and is NEVER checked!
     */
    static final int[] ALL_LOGLEVELS =  
        new int[]
            {
                EXCEPTION,
                INIT,
                PROGRESS,
                TIME_STEP,
                OUTPUT_STEP,
                METHOD_CALL,
                V_SCALE,
                QUANTITY,
                PARALLEL,
                A_IJ,
                PARTICLE,
                SCMVV
            };
    /**
     * Log level count
     */
    static final int LOGLEVEL_COUNT = ALL_LOGLEVELS.length;
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
    boolean isLogLevel(int aLogLevel);
    
    /**
     * Starts logger
     */
    void start();

    // <editor-fold defaultstate="collapsed" desc="LogLevel EXCEPTION related methods">
    /**
     * Appends with log level EXCEPTION
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMethod  Method in which exception occurred
     * @param aStacktrace Stack trace of excpetion
     */
    public void appendException(String aMethod, String aStacktrace);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel INIT related methods">
    /**
     * Appends with log level INIT
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationInit(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel PROGRESS related methods">
    /**
     * Appends with log level PROGRESS
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationProgress(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel TIME_STEP related methods">
    /**
     * Appends with log level TIME_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendTimeStep(String aMessage);

    /**
     * Appends start of time step with log level TIME_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepStart(int aTimeStep, long anId);

    /**
     * Appends end of time step with log level TIME_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepEnd(int aTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel OUTPUT_STEP related methods">
    /**
     * Appends with log level OUTPUT_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendOutputTimeStep(String aMessage);

    /**
     * Appends start of output time step with log level OUTPUT_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepStart(int anOutputTimeStep, long anId);

    /**
     * Appends end of output time step with log level OUTPUT_STEP
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepEnd(int anOutputTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel METHOD_CALL related methods">
    /**
     * Appends with log level METHOD_CALL
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendMethodCall(String aMessage);

    /**
     * Appends start of method with log level METHOD_CALL
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallStart(String aMessage, long anId);

    /**
     * Appends end of method with log level METHOD_CALL
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallEnd(String aMessage, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel V_SCALE related methods">
    /**
     * Appends with log level V_SCALE
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendVelocityScaleFactor(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel QUANTITY related methods">
    /**
     * Appends with log level QUANTITY
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendIntermediateResults(String aMessage);
    /**
     * Appends with log level QUANTITY
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anArrayName Array name
     * @param anArray Array
     */
    void appendIntermediateResults(String anArrayName, double[] anArray);    
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel PARALLEL related methods">
    /**
     * Appends with log level PARALLEL
 NOTE: Use isLogLevel() to check adequacy of append operation.
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
    // <editor-fold defaultstate="collapsed" desc="LogLevel PARTICLE related methods">
    /**
     * Appends with log level PARTICLE
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendParticleForceMagnitude(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel SCMVV related methods">
    /**
     * Appends with log level PARTICLE
 NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendScmvv(String aMessage);

    /**
     * True: SCMVV information accumulation, false: Otherwise
      * 
     * @param aValue Value
     */
    void setScmvvInformationAccumulation(boolean aValue);
    
    /**
     * True: SCMVV information accumulation, false: Otherwise
      * 
     * @return True: SCMVV information accumulation, false: Otherwise
     */
    boolean isScmvvInformationAccumulation();
    
    /**
     * SCMVV dissipative force particle index pair counter
     * 
     * @return SCMVV dissipative force particle index pair counter
     */
    AtomicInteger getScmvvFdissParticleIndexPairCounter();
    
    /**
     * SCMVV dissipative force particle index pair counter
     * 
     * @param aParticleIndexPairCounter Particle index pair counter
     */
    void setScmvvFdissParticleIndexPairCounter(AtomicInteger aParticleIndexPairCounter);

    /**
     * SCMVV dissipative force Rij_x adder
     * 
     * @return SCMVV dissipative force Rij_x adder
     */
    DoubleAdder getScmvvFdissRij_x_Adder();
    
    /**
     * SCMVV dissipative force Rij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissRij_x_Adder(DoubleAdder aDoubleAdder);
    
    /**
     * SCMVV dissipative force Abs(Rij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Rij_x) adder
     */
    DoubleAdder getScmvvFdissAbsRij_x_Adder();
    
    /**
     * SCMVV dissipative force Abs(Rij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissAbsRij_x_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force Vij_x adder
     * 
     * @return SCMVV dissipative force Vij_x adder
     */
    DoubleAdder getScmvvFdissVij_x_Adder();
    
    /**
     * SCMVV dissipative force Vij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissVij_x_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force Abs(Vij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Vij_x) adder
     */
    DoubleAdder getScmvvFdissAbsVij_x_Adder();
    
    /**
     * SCMVV dissipative force Abs(Vij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissAbsVij_x_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force Rij_x*Vij_x adder
     * 
     * @return SCMVV dissipative force Rij_x*Vij_x adder
     */
    DoubleAdder getScmvvFdissRijVij_x_Adder();
    
    /**
     * SCMVV dissipative force Rij_x*Vij_x adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissRijVij_x_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force Abs(Rij_x*Vij_x) adder
     * 
     * @return SCMVV dissipative force Abs(Rij_x*Vij)_x adder
     */
    DoubleAdder getScmvvFdissAbsRijVij_x_Adder();
    
    /**
     * SCMVV dissipative force Abs(Rij_x*Vij_x) adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissAbsRijVij_x_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force VdotR adder
     * 
     * @return SCMVV dissipative force VdotR adder
     */
    DoubleAdder getScmvvFdissVdotR_Adder();
    
    /**
     * SCMVV dissipative force VdotR adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissVdotR_Adder(DoubleAdder aDoubleAdder);

    /**
     * SCMVV dissipative force GammaFactor adder
     * 
     * @return SCMVV dissipative force GammaFactor adder
     */
    DoubleAdder getScmvvFdissGammaFactor_Adder();
    
    /**
     * SCMVV dissipative force GammaFactor adder
     * 
     * @param aDoubleAdder Double adder
     */
    void setScmvvFdissGammaFactor_Adder(DoubleAdder aDoubleAdder);
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
