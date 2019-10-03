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

/**
 * Group of unsafe (NON thread-safe) double adders
 * 
 * @author Achim Zielesny
 */
public class AdderGroup {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Potential energy double adder
     */
    private final DoubleAdderUnsafe potentialEnergyAdder;

    /**
     * Pressure tensor diagonal x term double adder
     */
    private final DoubleAdderUnsafe pressureXAdder;

    /**
     * Pressure tensor diagonal y term double adder
     */
    private final DoubleAdderUnsafe pressureYAdder;

    /**
     * Pressure tensor diagonal z term double adder
     */
    private final DoubleAdderUnsafe pressureZAdder;
    
    /**
     * Adder for G value of PNHLN integration type
     */
    private final DoubleAdderUnsafe pnhln_G_Adder;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public AdderGroup() {
        this.potentialEnergyAdder = new DoubleAdderUnsafe();
        this.pressureXAdder = new DoubleAdderUnsafe();
        this.pressureYAdder = new DoubleAdderUnsafe();
        this.pressureZAdder = new DoubleAdderUnsafe();
        this.pnhln_G_Adder = new DoubleAdderUnsafe();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Resets adders
     */
    public void reset() {
        this.potentialEnergyAdder.reset();
        this.pressureXAdder.reset();
        this.pressureYAdder.reset();
        this.pressureZAdder.reset();
        this.pnhln_G_Adder.reset();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Potential energy double adder
     * 
     * @return Potential energy double adder
     */
    public DoubleAdderUnsafe getPotentialEnergyAdder() {
        return this.potentialEnergyAdder;
    }

    /**
     * Pressure tensor diagonal x term double adder
     * 
     * @return Pressure tensor diagonal x term double adder
     */
    public DoubleAdderUnsafe getPressureXAdder() {
        return this.pressureXAdder;
    }

    /**
     * Pressure tensor diagonal y term double adder
     * 
     * @return Pressure tensor diagonal y term double adder
     */
    public DoubleAdderUnsafe getPressureYAdder() {
        return this.pressureYAdder;
    }

    /**
     * Pressure tensor diagonal z term double adder
     * 
     * @return Pressure tensor diagonal z term double adder
     */
    public DoubleAdderUnsafe getPressureZAdder() {
        return this.pressureZAdder;
    }

    /**
     * Adder for G value of PNHLN integration type
     * 
     * @return Adder for G value of PNHLN integration type
     */
    public DoubleAdderUnsafe getPnhln_G_Adder() {
        return this.pnhln_G_Adder;
    }
    // </editor-fold>

}
