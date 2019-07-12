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
 * Molecule acceleration description
 * 
 * @author Achim Zielesny
 */
public class MoleculeAccelerationDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;

    /**
     * x-component of acceleration
     */
    private final double accelerationX;

    /**
     * y-component of acceleration
     */
    private final double accelerationY;

    /**
     * z-component of acceleration
     */
    private final double accelerationZ;
    
    /**
     * Frequency
     */
    private final int frequency;
    
    /**
     * Maximum time step where acceleration is applied
     */
    private final int maxTimeStep;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * Note: NO checks are performed.
     * 
     * @param aMoleculeName Name of molecule
     * @param anAccelerationX x-component of acceleration
     * @param anAccelerationY y-component of acceleration
     * @param anAccelerationZ z-component of acceleration
     * @param aFrequency Frequency
     * @param aMaxTimeStep Maximum time step where acceleration is applied
     */
    public MoleculeAccelerationDescription(
        String aMoleculeName,
        double anAccelerationX,
        double anAccelerationY,
        double anAccelerationZ,
        int aFrequency,
        int aMaxTimeStep
    ) {
        this.moleculeName = aMoleculeName;
        this.accelerationX = anAccelerationX;
        this.accelerationY = anAccelerationY;
        this.accelerationZ = anAccelerationZ;
        this.frequency = aFrequency;
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
     * x-component of acceleration
     * 
     * @return x-component of acceleration
     */
    public final double getAccelerationX() {
        return this.accelerationX;
    }

    /**
     * y-component of acceleration
     * 
     * @return y-component of acceleration
     */
    public final double getAccelerationY() {
        return this.accelerationY;
    }

    /**
     * z-component of acceleration
     * 
     * @return z-component of acceleration
     */
    public final double getAccelerationZ() {
        return this.accelerationZ;
    }

    /**
     * Frequency
     * 
     * @return Frequency
     */
    public final int getFrequency() {
        return this.frequency;
    }

    /**
     * Maximum time step where acceleration is applied
     * 
     * @return Maximum time step where acceleration is applied
     */
    public final int getMaxTimeStep() {
        return this.maxTimeStep;
    }
    // </editor-fold>
    
}
