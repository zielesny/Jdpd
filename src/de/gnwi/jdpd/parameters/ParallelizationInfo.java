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
package de.gnwi.jdpd.parameters;

/**
 * Parallelisation info
 * 
 * @author Achim Zielesny
 */
public class ParallelizationInfo {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Minimum number of cells for parallelised task
     */
    private final int minimumParallelTaskCellNumber;

    /**
     * Minimum number of harmonic bonds for parallelised task
     */
    private final int minimumParallelTaskHarmonicBondNumber;
    
    /**
     * Number of parallel tasks for calculations
     */
    private final int parallelTaskNumber;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Maximum used parallel task number
     */
    private int maximumUsedParallelTaskNumber;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aMinimumParallelTaskCellNumber Minimum number of cells for parallelised task
     * @param aMinimumParallelTaskHarmonicBondNumber Minimum number of harmonic bonds for parallelised task
     * @param aParallelTaskNumber Number of parallel tasks for calculations
     * @throws IllegalArgumentException Thrown if argument is invalid
     */
    public ParallelizationInfo(
        int aMinimumParallelTaskCellNumber,
        int aMinimumParallelTaskHarmonicBondNumber,
        int aParallelTaskNumber) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aMinimumParallelTaskCellNumber < 2) {
            throw new IllegalArgumentException("ParallelizationInfo.Constructor: aMinimumParallelTaskCellNumber is less than 2.");
        }
        if (aMinimumParallelTaskHarmonicBondNumber < 2) {
            throw new IllegalArgumentException("ParallelizationInfo.Constructor: aMinimumParallelTaskHarmonicBondNumber is less than 2.");
        }
        if (aParallelTaskNumber < 1) {
            throw new IllegalArgumentException("ParallelizationInfo.Constructor: aParallelTaskNumber is less than 1.");
        }
        // </editor-fold>
        this.minimumParallelTaskCellNumber = aMinimumParallelTaskCellNumber;
        this.minimumParallelTaskHarmonicBondNumber = aMinimumParallelTaskHarmonicBondNumber;
        this.parallelTaskNumber = aParallelTaskNumber;
        
        this.maximumUsedParallelTaskNumber = 0;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Maximum used parallel task number
     * 
     * @param aValue Value
     */
    public void setMaximumUsedParallelTaskNumber(int aValue) {
        if (aValue > this.maximumUsedParallelTaskNumber) {
            this.maximumUsedParallelTaskNumber = aValue;
        }
    }

    /**
     * Maximum used parallel task number
     * 
     * @return Maximum used parallel task number
     */
    public int getMaximumUsedParallelTaskNumber() {
            return this.maximumUsedParallelTaskNumber;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Minimum number of cells for parallelised task
     * 
     * @return Minimum number of cells for parallelised task
     */
    public int getMinimumParallelTaskCellNumber() {
        return this.minimumParallelTaskCellNumber;
    }

    /**
     * Minimum number of harmonic bonds for parallelised task
     * 
     * @return Minimum number of harmonic bonds for parallelised task
     */
    public int getMinimumParallelTaskHarmonicBondNumber() {
        return this.minimumParallelTaskHarmonicBondNumber;
    }

    /**
     * Number of parallel tasks for calculations
     * 
     * @return Number of parallel tasks for calculations
     */
    public int getParallelTaskNumber() {
        return this.parallelTaskNumber;
    }
    // </editor-fold>
    
}
