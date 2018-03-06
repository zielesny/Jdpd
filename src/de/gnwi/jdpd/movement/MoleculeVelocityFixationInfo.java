/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2018  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.movement;

/**
 * Molecule velocity fixation info
 * 
 * @author Achim Zielesny
 */
public class MoleculeVelocityFixationInfo extends MoleculeVelocityFixationDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * First index in particle arrays for molecule
     */
    private final int firstIndex;

    /**
     * Exclusive last index in particle arrays for molecule
     */
    private final int exclusiveLastIndex;
    
    /**
     * Length in particle arrays for molecule
     */
    private final int length;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * Note: NO checks are performed.
     * 
     * @param aMoleculeName Name of molecule
     * @param anIsFixedX True: x-component of all particles of molecule is fixed, false: Otherwise
     * @param anIsFixedY True: y-component of all particles of molecule is fixed, false: Otherwise
     * @param anIsFixedZ True: z-component of all particles of molecule is fixed, false: Otherwise
     * @param aVelocityX x-component of velocity
     * @param aVelocityY y-component of velocity
     * @param aVelocityZ z-component of velocity
     * @param aFirstIndex First index in particle arrays for molecule
     * @param aLastIndex Last index in particle arrays for molecule
     */
    public MoleculeVelocityFixationInfo(
        String aMoleculeName,
        boolean anIsFixedX,
        boolean anIsFixedY,
        boolean anIsFixedZ,
        double aVelocityX,
        double aVelocityY,
        double aVelocityZ,
        int aFirstIndex,
        int aLastIndex) {
        super(
            aMoleculeName,
            anIsFixedX,
            anIsFixedY,
            anIsFixedZ,
            aVelocityX,
            aVelocityY,
            aVelocityZ
        );
        this.firstIndex = aFirstIndex;
        this.exclusiveLastIndex = aLastIndex + 1;
        // aLastIndex = 5, aFirstIndex = 3 -> length = 5 - 3 + 1 = 3
        this.length = aLastIndex - aFirstIndex + 1;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * First index in particle arrays for molecule
     * 
     * @return First index in particle arrays for molecule
     */
    public int getFirstIndex() {
        return this.firstIndex;
    }

    /**
     * Exclusive last index in particle arrays for molecule
     * 
     * @return Exclusive last index in particle arrays for molecule
     */
    public int getExclusiveLastIndex() {
        return this.exclusiveLastIndex;
    }
    
    /**
     * Length in particle arrays for molecule
     * 
     * @return Length in particle arrays for molecule
     */
    public int getLength() {
        return this.length;
    }
    // </editor-fold>
    
}
