/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpdsp.rg;

/**
 * Molecule particle radius of gyration (Rg) calculation description
 * 
 * @author Achim Zielesny
 */
public class MoleculeRgCalculationDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aMoleculeName Name of molecule
     */
    public MoleculeRgCalculationDescription(
        String aMoleculeName) {
        this.moleculeName = aMoleculeName;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Name of molecule
     * 
     * @return Name of molecule
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }
    // </editor-fold>
    
}
