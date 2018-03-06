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
package de.gnwi.jdpd.utilities;

/**
 * Gravitational acceleration
 * 
 * @author Achim Zielesny
 */
public class GravitationalAcceleration {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Gravitational acceleration in x-direction
     */
    private final double gravitationalAccelerationX;

    /**
     * Gravitational acceleration in y-direction
     */
    private final double gravitationalAccelerationY;

    /**
     * Gravitational acceleration in z-direction
     */
    private final double gravitationalAccelerationZ;
    
    /**
     * True: Gravitational acceleration, false: Otherwise
     */
    private final boolean isGravitationalAcceleration;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aGravitationalAccelerationX Gravitational acceleration in x-direction
     * @param aGravitationalAccelerationY Gravitational acceleration in y-direction
     * @param aGravitationalAccelerationZ Gravitational acceleration in z-direction
     */
    public GravitationalAcceleration(double aGravitationalAccelerationX, double aGravitationalAccelerationY, double aGravitationalAccelerationZ) {
        this.gravitationalAccelerationX = aGravitationalAccelerationX;
        this.gravitationalAccelerationY = aGravitationalAccelerationY;
        this.gravitationalAccelerationZ = aGravitationalAccelerationZ;
        this.isGravitationalAcceleration = this.gravitationalAccelerationX != 0.0 || this.gravitationalAccelerationY != 0.0 || this.gravitationalAccelerationZ != 0.0; 
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Gravitational acceleration in x-direction
     * 
     * @return Gravitational acceleration in x-direction
     */
    public double getGravitationalAccelerationX() {
        return this.gravitationalAccelerationX;
    }

    /**
     * Gravitational acceleration in y-direction
     * 
     * @return Gravitational acceleration in y-direction
     */
    public double getGravitationalAccelerationY() {
        return this.gravitationalAccelerationY;
    }

    /**
     * Gravitational acceleration in z-direction
     * 
     * @return Gravitational acceleration in z-direction
     */
    public double getGravitationalAccelerationZ() {
        return this.gravitationalAccelerationZ;
    }
    
    /**
     * True: Gravitational acceleration, false: Otherwise
     * 
     * @return True: Gravitational acceleration, false: Otherwise
     */
    public boolean isGravitationalAcceleration() {
        return this.isGravitationalAcceleration;
    }
    // </editor-fold>
    
}

