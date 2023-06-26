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
 * Interface for calculators
 * 
 * @author Achim Zielesny
 */
public interface ICalculator {

    // <editor-fold defaultstate="collapsed" desc="Methods">
    /**
     * Executor service shutdown
     */
    void shutdownExecutorService();
    
    /**
     * Return accumulated (total) sum of all potential energy adders
     * 
     * @return Return accumulated (total) sum of all potential energy adders
     */
    double getAccumulatedPotentialEnergyAddersSum();
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal x term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal x term adders
     */
    double getAccumulatedPressureXAddersSum();
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal y term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal y term adders
     */
    double getAccumulatedPressureYAddersSum();
    
    /**
     * Return accumulated (total) sum of all pressure tensor diagonal z term adders
     * (No checks are performed)
     * 
     * @return Return accumulated (total) sum of all pressure tensor diagonal z term adders
     */
    double getAccumulatedPressureZAddersSum();
    // </editor-fold>
    
}
