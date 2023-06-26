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
package de.gnwi.jdpd.samples;

import de.gnwi.jdpd.interfaces.IProgressMonitor;

/**
 * Progress monitor
 * 
 * @author Achim Zielesny
 */
public class ProgressMonitor implements IProgressMonitor {

    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Simulation state
     */
    private SimulationState simulationState;
    
    /**
     * Progress in percent
     */
    private int progressInPercent;

    /**
     * Remaining time info
     */
    private String remainingTimeInfo;
    
    /**
     * True: simulation was stopped, false: Otherwise
     */
    private boolean wasSimulationStopped;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public ProgressMonitor() {
        this.simulationState = SimulationState.NOT_STARTED;
        this.progressInPercent = -1;
        this.remainingTimeInfo = null;
        this.wasSimulationStopped = false;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties">
    /**
     * Sets simulation state information
     * 
     * @param aSimulationState Simulation state
     */
    @Override
    public void setSimulationState(SimulationState aSimulationState) {
        this.simulationState = aSimulationState;
    }
    
    /**
     * Returns simulation state information
     * 
     * @return Simulation state information
     */
    @Override
    public SimulationState getSimulationState() {
        return this.simulationState;
    }
    
    /**
     * Sets progress of simulation during simulation state TIME_STEP_INTEGRATION
     * in percent
     * 
     * @param aProgressInPercent Progress in percent
     */
    @Override
    public void setProgressInPercent(int aProgressInPercent) {
        this.progressInPercent = aProgressInPercent;
    }
    
    /**
     * Returns progress of simulation during simulation state 
     * TIME_STEP_INTEGRATION in percent
     * 
     * @return Progress of simulation during simulation state 
     * TIME_STEP_INTEGRATION in percent
     */
    @Override
    public int getProgressInPercent() {
        return this.progressInPercent;
    }
    
    /**
     * Sets string with remaining time information
     * 
     * @param aRemainingTimeInfo Remaining time info
     */
    @Override
    public void setRemainingTime(String aRemainingTimeInfo) {
        this.remainingTimeInfo = aRemainingTimeInfo;
    }
    
    /**
     * Returns string with remaining time information
     * 
     * @return String with remaining time information
     */
    @Override
    public String getRemainingTime() {
        return this.remainingTimeInfo;
    }
    
    /**
     * Returns if simulation has finished.
     * 
     * @return True: Simulation has finished, false: Otherwise
     */
    @Override
    public boolean hasFinished() {
        return this.simulationState == SimulationState.FINISHED_WITH_FAILURE
            || this.simulationState == SimulationState.FINISHED_WITH_SUCCESS;
    }
    
    /**
     * Sets flag that simulation was stopped.
     */
    @Override
    public void setSimulationStoppedFlag() {
        this.wasSimulationStopped = true;
    }
    
    /**
     * Returns if simulation was interrupted.
     * 
     * @return True: Simulation was interrupted, false: Otherwise
     */
    @Override
    public boolean wasInterrupted() {
        return this.wasSimulationStopped;
    }
    // </editor-fold>
    
}
