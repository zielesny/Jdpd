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
package de.gnwi.jdpd.movement;

/**
 * Molecule sphere info
 * 
 * @author Achim Zielesny
 */
public class MoleculeSphereInfo extends MoleculeSphereDescription {
    
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
     * (No checks are performed)
     * 
     * @param aMoleculeName Name of molecule
     * @param anIsExclusiveSphere True: Molecule is not allowed inside volume of 
     * sphere and must be outside of the sphere , false: Molecule is not allowed 
     * outside sphere and must be inside the sphere
     * @param aSphereCenterX X component of sphere center point
     * @param aSphereCenterY Y component of sphere center point
     * @param aSphereCenterZ Z component of sphere center point
     * @param aSphereRadius Radius of sphere
     * @param aMaxTimeStep Maximum time step for application
     * @param aFirstIndex First index in particle arrays for molecule
     * @param aLastIndex Last index in particle arrays for molecule
     */
    public MoleculeSphereInfo(
        String aMoleculeName,
        boolean anIsExclusiveSphere,
        double aSphereCenterX,
        double aSphereCenterY,
        double aSphereCenterZ,
        double aSphereRadius,
        int aMaxTimeStep,
        int aFirstIndex,
        int aLastIndex
    ) {
        super(
            aMoleculeName,
            anIsExclusiveSphere,
            aSphereCenterX,
            aSphereCenterY,
            aSphereCenterZ,
            aSphereRadius,
            aMaxTimeStep
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
