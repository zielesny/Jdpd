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
package de.gnwi.jdpdsp.interfaces;

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Interface for simulation logger
 * 
 * @author Achim Zielesny
 */
public interface ILogger {

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
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMethod  Method in which exception occurred
     * @param aStacktrace Stack trace of excpetion
     */
    public void appendException(String aMethod, String aStacktrace);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel INIT related methods">
    /**
     * Appends with log level INIT
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationInit(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel PROGRESS related methods">
    /**
     * Appends with log level PROGRESS
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendSimulationProgress(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel TIME_STEP related methods">
    /**
     * Appends with log level TIME_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendTimeStep(String aMessage);

    /**
     * Appends start of time step with log level TIME_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepStart(int aTimeStep, long anId);

    /**
     * Appends end of time step with log level TIME_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aTimeStep Time step
     * @param anId ID for method identification
     */
    void appendTimeStepEnd(int aTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel OUTPUT_STEP related methods">
    /**
     * Appends with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendOutputTimeStep(String aMessage);

    /**
     * Appends start of output time step with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepStart(int anOutputTimeStep, long anId);

    /**
     * Appends end of output time step with log level OUTPUT_STEP
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anOutputTimeStep Output time step
     * @param anId ID for method identification
     */
    void appendOutputTimeStepEnd(int anOutputTimeStep, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel METHOD_CALL related methods">
    /**
     * Appends with log level METHOD_CALL
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendMethodCall(String aMessage);

    /**
     * Appends start of method with log level METHOD_CALL
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallStart(String aMessage, long anId);

    /**
     * Appends end of method with log level METHOD_CALL
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     * @param anId ID for method identification
     */
    void appendMethodCallEnd(String aMessage, long anId);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel V_SCALE related methods">
    /**
     * Appends with log level V_SCALE
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendVelocityScaleFactor(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel QUANTITY related methods">
    /**
     * Appends with log level QUANTITY
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendIntermediateResults(String aMessage);
    /**
     * Appends with log level QUANTITY
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param anArrayName Array name
     * @param anArray Array
     */
    void appendIntermediateResults(String anArrayName, float[] anArray);    
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel PARALLEL related methods">
    /**
     * Appends with log level PARALLEL
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
    // <editor-fold defaultstate="collapsed" desc="LogLevel PARTICLE related methods">
    /**
     * Appends with log level PARTICLE
     * NOTE: Use isLogLevel() to check adequacy of append operation.
     * 
     * @param aMessage Message
     */
    void appendParticleForceMagnitude(String aMessage);
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="LogLevel SCMVV related methods">
    /**
     * Appends with log level PARTICLE
     * NOTE: Use isLogLevel() to check adequacy of append operation.
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
