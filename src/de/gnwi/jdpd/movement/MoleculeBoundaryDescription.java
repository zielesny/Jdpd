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
 * Molecule boundary description
 * 
 * @author Achim Zielesny
 */
public class MoleculeBoundaryDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;

    /**
     * True: Boundary in x-direction is defined, false: Otherwise
     */
    private final boolean isBoundaryX;
    
    /**
     * Boundary x-min in simulation box
     */
    private final double xMin;

    /**
     * Boundary x-max in simulation box
     */
    private final double xMax;

    /**
     * True: Boundary in y-direction is defined, false: Otherwise
     */
    private final boolean isBoundaryY;

    /**
     * Boundary y-min in simulation box
     */
    private final double yMin;

    /**
     * Boundary y-max in simulation box
     */
    private final double yMax;

    /**
     * True: Boundary in z-direction is defined, false: Otherwise
     */
    private final boolean isBoundaryZ;
    
    /**
     * Boundary z-min in simulation box
     */
    private final double zMin;

    /**
     * Boundary z-max in simulation box
     */
    private final double zMax;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * Note: NO checks are performed.
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
     */
    public MoleculeBoundaryDescription(
        String aMoleculeName,
        boolean anIsBoundaryX,
        double aXmin,
        double aXmax,
        boolean anIsBoundaryY,
        double aYmin,
        double aYmax,
        boolean anIsBoundaryZ,
        double aZmin,
        double aZmax
    ) {
        this.moleculeName = aMoleculeName;
        this.isBoundaryX = anIsBoundaryX;
        this.xMin = aXmin;
        this.xMax = aXmax;
        this.isBoundaryY = anIsBoundaryY;
        this.yMin = aYmin;
        this.yMax = aYmax;
        this.isBoundaryZ = anIsBoundaryZ;
        this.zMin = aZmin;
        this.zMax = aZmax;
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

    /**
     * True: Boundary in x-direction is defined, false: Otherwise
     * 
     * @return True: Boundary in x-direction is defined, false: Otherwise
     */
    public boolean isBoundaryX() {
        return isBoundaryX;
    }
    
    /**
     * Boundary x-min in simulation box
     * 
     * @return Boundary x-min in simulation box
     */
    public double getXmin() {
        return this.xMin;
    }

    /**
     * Boundary x-max in simulation box
     * 
     * @return Boundary x-max in simulation box
     */
    public double getXmax() {
        return this.xMax;
    }

    /**
     * True: Boundary in y-direction is defined, false: Otherwise
     * 
     * @return True: Boundary in y-direction is defined, false: Otherwise
     */
    public boolean isBoundaryY() {
        return isBoundaryY;
    }

    /**
     * Boundary y-min in simulation box
     * 
     * @return Boundary y-min in simulation box
     */
    public double getYmin() {
        return this.yMin;
    }

    /**
     * Boundary y-max in simulation box
     * 
     * @return Boundary y-max in simulation box
     */
    public double getYmax() {
        return this.yMax;
    }

    /**
     * True: Boundary in z-direction is defined, false: Otherwise
     * 
     * @return True: Boundary in z-direction is defined, false: Otherwise
     */
    public boolean isBoundaryZ() {
        return isBoundaryZ;
    }

    /**
     * Boundary z-min in simulation box
     * 
     * @return Boundary z-min in simulation box
     */
    public double getZmin() {
        return this.zMin;
    }

    /**
     * Boundary z-max in simulation box
     * 
     * @return Boundary z-max in simulation box
     */
    public double getZmax() {
        return this.zMax;
    }
    // </editor-fold>
    
}
