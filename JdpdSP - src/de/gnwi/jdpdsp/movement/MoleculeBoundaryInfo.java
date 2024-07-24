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
package de.gnwi.jdpdsp.movement;

/**
 * Molecule boundary info
 * 
 * @author Achim Zielesny
 */
public class MoleculeBoundaryInfo extends MoleculeBoundaryDescription {
    
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
     * @param anIsBoundaryX True: Boundary in x-direction is defined, false: Otherwise
     * @param aXmin Boundary x-min in simulation box
     * @param aXmax Boundary x-max in simulation box
     * @param anIsBoundaryY True: Boundary in y-direction is defined, false: Otherwise
     * @param aYmin Boundary y-min in simulation box
     * @param aYmax Boundary y-max in simulation box
     * @param anIsBoundaryZ True: Boundary in z-direction is defined, false: Otherwise
     * @param aZmin Boundary z-min in simulation box
     * @param aZmax Boundary z-max in simulation box
     * @param aMaxTimeStep Maximum time step for application
     * @param aFirstIndex First index in particle arrays for molecule
     * @param aLastIndex Last index in particle arrays for molecule
     */
    public MoleculeBoundaryInfo(
        String aMoleculeName,
        boolean anIsBoundaryX,
        float aXmin,
        float aXmax,
        boolean anIsBoundaryY,
        float aYmin,
        float aYmax,
        boolean anIsBoundaryZ,
        float aZmin,
        float aZmax,
        int aMaxTimeStep,
        int aFirstIndex,
        int aLastIndex
    ) {
        super(
            aMoleculeName,
            anIsBoundaryX,
            aXmin,
            aXmax,
            anIsBoundaryY,
            aYmin,
            aYmax,
            anIsBoundaryZ,
            aZmin,
            aZmax,
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
