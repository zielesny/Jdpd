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
package de.gnwi.jdpdsp.utilities;

import java.util.concurrent.ExecutorService;

/**
 * Utility class for calculators
 * 
 * @author Achim Zielesny
 */
public class CalculatorUtils {

    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Executor service
     */
    private ExecutorService executorService;
    
    /**
     * Adder groups
     */
    private AdderGroup[] adderGroups;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public CalculatorUtils() {
        this.executorService = null;
        this.adderGroups = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Resets adders
     * (No checks are performed)
     */
    public void resetAdders() {
        if (this.adderGroups != null) {
            for (AdderGroup tmpAdderGroup : this.adderGroups) {
                tmpAdderGroup.reset();
            }
        }
    }

    /**
     * Executor service shutdown
     */
    public void shutdownExecutorService() {
        if (this.executorService != null) {
            this.executorService.shutdown();
            this.executorService = null;
        }
    }
    
    /**
     * Return accumulated (total) sum of all potential energy adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all potential energy adders
     */
    public float getAccumulatedPotentialEnergyAddersSum() {
        float tmpTotalSum = 0.0f;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPotentialEnergyAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal x term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal x term adders
     */
    public float getAccumulatedPressureXAddersSum() {
        float tmpTotalSum = 0.0f;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureXAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal y term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal y term adders
     */
    public float getAccumulatedPressureYAddersSum() {
        float tmpTotalSum = 0.0f;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureYAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal z term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal z term adders
     */
    public float getAccumulatedPressureZAddersSum() {
        float tmpTotalSum = 0.0f;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPressureZAdder().getSum();
        }
        return tmpTotalSum;
    }
    
    /**
     * Return accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all adders for G value of PNHLN 
     * integration
     */
    public float getAccumulatedPnhln_G_AddersSum() {
        float tmpTotalSum = 0.0f;
        for (AdderGroup tmpAdderGroup : this.adderGroups) {
            tmpTotalSum += tmpAdderGroup.getPnhln_G_Adder().getSum();
        }
        return tmpTotalSum;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get/set)">
    /**
     * Executor service
     * 
     * @return Executor service
     */
    public ExecutorService getExecutorService() {
        return this.executorService;
    }

    /**
     * Executor service
     * (No checks are performed)
     * 
     * @param anExecutorService Executer service
     */
    public void setExecutorService(ExecutorService anExecutorService) {
        this.executorService =  anExecutorService;
    }
    
    /**
     * Adder groups
     * 
     * @return Adder groups
     */
    public AdderGroup[] getAdderGroups() {
        return this.adderGroups;
    }
    
    /**
     * Adder groups
     * (No checks are performed)
     * 
     * @param anAdderGroups Adder groups
     */
    public void setAdderGroups(AdderGroup[] anAdderGroups) {
        this.adderGroups = anAdderGroups;
    }
    // </editor-fold>
    
}
