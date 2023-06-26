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
 * Molecule sphere description
 * 
 * @author Achim Zielesny
 */
public class MoleculeSphereDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;

    /**
     * True: Molecule is not allowed inside volume of sphere and must be outside
     * of the sphere , false: Molecule is not allowed outside sphere and must 
     * be inside the sphere
     */
    private final boolean isExclusiveSphere;
    
    /**
     * X component of sphere center point
     */
    private final double sphereCenterX;

    /**
     * Y component of sphere center point
     */
    private final double sphereCenterY;

    /**
     * Z component of sphere center point
     */
    private final double sphereCenterZ;

    /**
     * Radius of sphere
     */
    private final double sphereRadius;

    /**
     * Diameter of sphere (= 2 * sphereRadius)
     */
    private final double sphereDiameter;

    /**
     * Square of radius of sphere
     */
    private final double sphereRadiusSquare;

    /**
     * Maximum time step for application
     */
    private final int maxTimeStep;
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
     */
    public MoleculeSphereDescription(
        String aMoleculeName,
        boolean anIsExclusiveSphere,
        double aSphereCenterX,
        double aSphereCenterY,
        double aSphereCenterZ,
        double aSphereRadius,
        int aMaxTimeStep
    ) {
        this.moleculeName = aMoleculeName;
        this.isExclusiveSphere = anIsExclusiveSphere;
        this.sphereCenterX = aSphereCenterX;
        this.sphereCenterY = aSphereCenterY;
        this.sphereCenterZ = aSphereCenterZ;
        this.sphereRadius = aSphereRadius;
        this.sphereDiameter = aSphereRadius + aSphereRadius;
        this.sphereRadiusSquare = aSphereRadius * aSphereRadius;
        this.maxTimeStep = aMaxTimeStep;
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
     * True: Molecule is not allowed inside volume of sphere and must be outside
     * of the sphere , false: Molecule is not allowed outside sphere and must 
     * be inside the sphere
     * 
     * @return True: Molecule is not allowed inside volume of sphere and must be 
     * outside of the sphere , false: Molecule is not allowed outside sphere and 
     * must be inside the sphere
     */
    public boolean isExclusiveSphere() {
        return this.isExclusiveSphere;
    }
    
    /**
     * X component of sphere center point
     * 
     * @return X component of sphere center point
     */
    public double getSphereCenterX() {
        return this.sphereCenterX;
    }
    
    /**
     * Y component of sphere center point
     * 
     * @return Y component of sphere center point
     */
    public double getSphereCenterY() {
        return this.sphereCenterY;
    }

    /**
     * Z component of sphere center point
     * 
     * @return Z component of sphere center point
     */
    public double getSphereCenterZ() {
        return this.sphereCenterZ;
    }
        
    /**
     * Radius of sphere
     * 
     * @return Radius of sphere
     */
    public double getSphereRadius() {
        return this.sphereRadius;
    }
            
    /**
     * Diameter of sphere
     * 
     * @return Diameter of sphere
     */
    public double getSphereDiameter() {
        return this.sphereDiameter;
    }

    /**
     * Square of radius of sphere
     * 
     * @return Square of radius of sphere
     */
    public double getSphereRadiusSquare() {
        return this.sphereRadiusSquare;
    }

    /**
     * Maximum time step for application
     * 
     * @return Maximum time step for application
     */
    public final int getMaxTimeStep() {
        return this.maxTimeStep;
    }
    // </editor-fold>
    
}