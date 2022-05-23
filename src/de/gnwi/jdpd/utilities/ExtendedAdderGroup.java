/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2022  Achim Zielesny (achim.zielesny@googlemail.com)
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
 * Extended group of Unsafe (NON thread-safe) double adders
 * 
 * @author Achim Zielesny
 */
public class ExtendedAdderGroup extends AdderGroup {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * DPD potential energy double adder
     */
    private final DoubleAdderUnsafe dpdPotentialEnergyAdder;

    /**
     * Bond potential energy double adder
     */
    private final DoubleAdderUnsafe bondPotentialEnergyAdder;

    /**
     * Electrostatics potential energy double adder
     */
    private final DoubleAdderUnsafe electrostaticsPotentialEnergyAdder;
    
    /**
     * DPD pressure tensor diagonal x term double adder
     */
    private final DoubleAdderUnsafe dpdPressureXAdder;

    /**
     * DPD pressure tensor diagonal y term double adder
     */
    private final DoubleAdderUnsafe dpdPressureYAdder;

    /**
     * DPD pressure tensor diagonal z term double adder
     */
    private final DoubleAdderUnsafe dpdPressureZAdder;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public ExtendedAdderGroup() {
        super();
        this.dpdPotentialEnergyAdder = new DoubleAdderUnsafe();
        this.bondPotentialEnergyAdder = new DoubleAdderUnsafe();
        this.electrostaticsPotentialEnergyAdder = new DoubleAdderUnsafe();
        this.dpdPressureXAdder = new DoubleAdderUnsafe();
        this.dpdPressureYAdder = new DoubleAdderUnsafe();
        this.dpdPressureZAdder = new DoubleAdderUnsafe();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Resets adders
     */
    @Override
    public void reset() {
        super.reset();
        this.dpdPotentialEnergyAdder.reset();
        this.bondPotentialEnergyAdder.reset();
        this.electrostaticsPotentialEnergyAdder.reset();
        this.dpdPressureXAdder.reset();
        this.dpdPressureYAdder.reset();
        this.dpdPressureZAdder.reset();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * DPD potential energy double adder
     * 
     * @return DPD potential energy double adder
     */
    public DoubleAdderUnsafe getDpdPotentialEnergyAdder() {
        return this.dpdPotentialEnergyAdder;
    }

    /**
     * Bond potential energy double adder
     * 
     * @return Bond potential energy double adder
     */
    public DoubleAdderUnsafe getBondPotentialEnergyAdder() {
        return this.bondPotentialEnergyAdder;
    }

    /**
     * Electrostatics potential energy double adder
     * 
     * @return Electrostatics potential energy double adder
     */
    public DoubleAdderUnsafe getElectrostaticsPotentialEnergyAdder() {
        return this.electrostaticsPotentialEnergyAdder;
    }

    /**
     * DPD pressure tensor diagonal x term double adder
     * 
     * @return DPD pressure tensor diagonal x term double adder
     */
    public DoubleAdderUnsafe getDpdPressureXAdder() {
        return this.dpdPressureXAdder;
    }

    /**
     * DPD pressure tensor diagonal y term double adder
     * 
     * @return DPD pressure tensor diagonal y term double adder
     */
    public DoubleAdderUnsafe getDpdPressureYAdder() {
        return this.dpdPressureYAdder;
    }

    /**
     * DPD pressure tensor diagonal z term double adder
     * 
     * @return DPD pressure tensor diagonal z term double adder
     */
    public DoubleAdderUnsafe getDpdPressureZAdder() {
        return this.dpdPressureZAdder;
    }
    // </editor-fold>

}
