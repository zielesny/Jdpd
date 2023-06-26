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

/**
 * Extended group of Unsafe (NON thread-safe) float adders
 * 
 * @author Achim Zielesny
 */
public class ExtendedAdderGroup extends AdderGroup {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * DPD potential energy float adder
     */
    private final FloatAdderUnsafe dpdPotentialEnergyAdder;

    /**
     * Bond potential energy float adder
     */
    private final FloatAdderUnsafe bondPotentialEnergyAdder;

    /**
     * Electrostatics potential energy float adder
     */
    private final FloatAdderUnsafe electrostaticsPotentialEnergyAdder;
    
    /**
     * DPD pressure tensor diagonal x term float adder
     */
    private final FloatAdderUnsafe dpdPressureXAdder;

    /**
     * DPD pressure tensor diagonal y term float adder
     */
    private final FloatAdderUnsafe dpdPressureYAdder;

    /**
     * DPD pressure tensor diagonal z term float adder
     */
    private final FloatAdderUnsafe dpdPressureZAdder;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public ExtendedAdderGroup() {
        super();
        this.dpdPotentialEnergyAdder = new FloatAdderUnsafe();
        this.bondPotentialEnergyAdder = new FloatAdderUnsafe();
        this.electrostaticsPotentialEnergyAdder = new FloatAdderUnsafe();
        this.dpdPressureXAdder = new FloatAdderUnsafe();
        this.dpdPressureYAdder = new FloatAdderUnsafe();
        this.dpdPressureZAdder = new FloatAdderUnsafe();
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
     * DPD potential energy float adder
     * 
     * @return DPD potential energy float adder
     */
    public FloatAdderUnsafe getDpdPotentialEnergyAdder() {
        return this.dpdPotentialEnergyAdder;
    }

    /**
     * Bond potential energy float adder
     * 
     * @return Bond potential energy float adder
     */
    public FloatAdderUnsafe getBondPotentialEnergyAdder() {
        return this.bondPotentialEnergyAdder;
    }

    /**
     * Electrostatics potential energy float adder
     * 
     * @return Electrostatics potential energy float adder
     */
    public FloatAdderUnsafe getElectrostaticsPotentialEnergyAdder() {
        return this.electrostaticsPotentialEnergyAdder;
    }

    /**
     * DPD pressure tensor diagonal x term float adder
     * 
     * @return DPD pressure tensor diagonal x term float adder
     */
    public FloatAdderUnsafe getDpdPressureXAdder() {
        return this.dpdPressureXAdder;
    }

    /**
     * DPD pressure tensor diagonal y term float adder
     * 
     * @return DPD pressure tensor diagonal y term float adder
     */
    public FloatAdderUnsafe getDpdPressureYAdder() {
        return this.dpdPressureYAdder;
    }

    /**
     * DPD pressure tensor diagonal z term float adder
     * 
     * @return DPD pressure tensor diagonal z term float adder
     */
    public FloatAdderUnsafe getDpdPressureZAdder() {
        return this.dpdPressureZAdder;
    }
    // </editor-fold>

}
