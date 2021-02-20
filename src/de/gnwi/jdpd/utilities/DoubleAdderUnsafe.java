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
package de.gnwi.jdpd.utilities;

/**
 * Double adder
 * NOTE: This implementation is NOT thread-safe.
 * 
 * @author Achim Zielesny
 */
public class DoubleAdderUnsafe {

    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Sum
     */
    private double sum;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public DoubleAdderUnsafe() {
        this.sum = 0.0;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Adds a value to the sum
     * 
     * @param aValue value
     */
    public void add(double aValue) {
        this.sum += aValue;
    }
    
    /**
     * Resets sum to 0.0
     */
    public void reset() {
        this.sum = 0.0;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Sum
     * 
     * @return Sum
     */
    public double getSum() {
        return this.sum;
    }
    // </editor-fold>

}
