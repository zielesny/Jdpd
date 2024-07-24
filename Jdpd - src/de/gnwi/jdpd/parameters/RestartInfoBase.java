/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpd.parameters;

/**
 * Restart information base
 * 
 * @author Achim Zielesny
 */
public class RestartInfoBase {
    
    // <editor-fold defaultstate="collapsed" desc="protected static final class variables">
    /**
     * Version
     */
    protected final static String VERSION_1_0_0 = "Version 1.0.0";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="protected final class variables">
    /**
     * Number of additional simulation time steps
     */
    protected final int additionalTimeStepNumber;
    
    /**
     * True: Velocity initialization is to be performed, false: Otherwise
     */
    protected final boolean isVelocityInitialization;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="protected class variables">
    /**
     * Last simulation time step to which the particle positions and velocities belong
     */
    protected int lastTimeStep;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aLastTimeStep Last simulation time step to which the particle 
     * positions and velocities belong
     * @param anAdditionalTimeStepNumber Number of additional simulation time 
     * @param anIsVelocityInitialization True: Velocity initialization is to be 
     * performed, false: Otherwise
     * steps
     */
    public RestartInfoBase(
        int aLastTimeStep,
        int anAdditionalTimeStepNumber,
        boolean anIsVelocityInitialization
    ) {
        this.lastTimeStep = aLastTimeStep;
        this.additionalTimeStepNumber = anAdditionalTimeStepNumber;
        this.isVelocityInitialization = anIsVelocityInitialization;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected properties (get)">
    /**
     * Last simulation time step to which the particle positions and velocities belong
     * 
     * @return Last simulation time step to which the particle positions and velocities belong
     */
    public int getLastTimeStep () {
        return this.lastTimeStep;
    }

    /**
     * Number of additional simulation time steps
     * 
     * @return Number of additional simulation time steps
     */
    public int getAdditionalTimeStepNumber() {
        return this.additionalTimeStepNumber;
    }

    /**
     * Flag for velocity initialization
     * 
     * @return True: Velocity initialization is to be performed, false: 
     * Otherwise
     */
    public boolean isVelocityInitialization() {
        return this.isVelocityInitialization;
    }
    // </editor-fold>
    
}
