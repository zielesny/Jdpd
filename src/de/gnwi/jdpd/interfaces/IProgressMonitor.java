/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2019  Achim Zielesny (achim.zielesny@googlemail.com)
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
 * Interface for progress monitor
 * 
 * @author Achim Zielesny
 */
public interface IProgressMonitor {

    // <editor-fold defaultstate="collapsed" desc="Enum SimulationState">
    /**
     * Simulation state
     */
    public enum SimulationState {
            
        /**
         * Simulation state NOT_STARTED
         */
        NOT_STARTED,
        /**
         * Simulation state STARTED
         */
        STARTED,
        /**
         * Simulation state PRE_PROCESSING
         */
        PRE_PROCESSING,
        /**
         * Simulation state TIME_STEP_INTEGRATION
         */
        TIME_STEP_INTEGRATION,
        /**
         * Simulation state POST_PROCESSING
         */
        POST_PROCESSING,
        /**
         * Simulation state FINISHED_WITH_SUCCESS
         */
        FINISHED_WITH_SUCCESS,
        /**
         * Simulation state FINISHED_WITH_FAILURE
         */
        FINISHED_WITH_FAILURE,

    }    
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public progress related properties">
    /**
     * Sets simulation state information
     * 
     * @param aSimulationState Simulation state
     */
    void setSimulationState(SimulationState aSimulationState);
    
    /**
     * Sets progress of simulation during simulation state TIME_STEP_INTEGRATION
     * in percent
     * 
     * @param aProgressInPercent Progress in percent
     */
    void setProgressInPercent(int aProgressInPercent); 
    
    /**
     * Sets string with remaining time information
     * 
     * @param aRemainingTimeInfo Remaining time info
     */
    void setRemainingTime(String aRemainingTimeInfo);
    
    /**
     * Returns if simulation has finished.
     * 
     * @return True: Simulation has finished, false: Otherwise
     */
    boolean hasFinished();

    /**
     * Sets flag that simulation was stopped.
     */
    void setSimulationStoppedFlag();
    
    /**
     * Returns if simulation was interrupted.
     * 
     * @return True: Simulation was interrupted, false: Otherwise
     */
    boolean wasInterrupted();
    // </editor-fold>
    
}
