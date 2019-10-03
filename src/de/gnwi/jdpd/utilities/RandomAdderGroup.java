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
package de.gnwi.jdpd.utilities;

import de.gnwi.jdpd.interfaces.IRandom;

/**
 * Group of unsafe (NON thread-safe) double adders plus random number generator
 * 
 * @author Achim Zielesny
 */
public class RandomAdderGroup extends AdderGroup {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Random number generator
     */
    private final IRandom randomNumberGenerator;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aRandomNumberGenerator Random number generator
     * @throws IllegalArgumentException Thrown if argument is illegal
     */
    public RandomAdderGroup(IRandom aRandomNumberGenerator) {
        super();
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aRandomNumberGenerator == null) {
            throw new IllegalArgumentException("UnsafeRandomAdders.Constructor: aRandomNumberGenerator is null.");
        }
        // </editor-fold>
        this.randomNumberGenerator = aRandomNumberGenerator;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Random number generator
     * 
     * @return Random number generator
     */
    public IRandom getRandomNumberGenerator() {
        return this.randomNumberGenerator;
    }
    // </editor-fold>

}
