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
package de.gnwi.jdpd.rg;

/**
 * Molecule particle radius of gyration (Rg) value
 * 
 * @author Achim Zielesny
 */
public class MoleculeRgValue {
    
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Molecule name
     */
    private String moleculeName;

    /**
     * Rg value
     */
    private double rgValue;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public MoleculeRgValue() {
    }

    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aMoleculeName Molecule name
     * @param aRgValue Rg value
     * 
     */
    public MoleculeRgValue(String aMoleculeName, double aRgValue) {
        this.setMoleculeRgValue(aMoleculeName, aRgValue);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Set molecule and Rg value
     * (No checks are performed)
     * 
     * @param aMoleculeName Molecule name
     * @param aRgValue Rg value
     */
    public final void setMoleculeRgValue(String aMoleculeName, double aRgValue) {
        this.moleculeName = aMoleculeName;
        this.rgValue = aRgValue;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Molecule name
     * 
     * @return Molecule name
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }

    /**
     * Rg value
     * 
     * @return Rg value
     */
    public double getRgValue() {
        return this.rgValue;
    }
    // </editor-fold>
}
