/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpd.interfaces;

/**
 * Interface for log levels
 * 
 * @author Achim Zielesny
 */
public interface ILogLevel {

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
    
}
