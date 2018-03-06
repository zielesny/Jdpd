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
 * Molecule velocity fixation description
 * 
 * @author Achim Zielesny
 */
public class MoleculeVelocityFixationDescription extends MoleculeFixationDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * x-component of velocity
     */
    private final double velocityX;

    /**
     * y-component of velocity
     */
    private final double velocityY;

    /**
     * z-component of velocity
     */
    private final double velocityZ;
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
     */
    public MoleculeVelocityFixationDescription(
        String aMoleculeName,
        boolean anIsFixedX,
        boolean anIsFixedY,
        boolean anIsFixedZ,
        double aVelocityX,
        double aVelocityY,
        double aVelocityZ) {
        super(
            aMoleculeName,
            anIsFixedX,
            anIsFixedY,
            anIsFixedZ
        );
        this.velocityX = aVelocityX;
        this.velocityY = aVelocityY;
        this.velocityZ = aVelocityZ;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * x-component of velocity
     * 
     * @return x-component of velocity
     */
    public final double getVelocityX() {
        return this.velocityX;
    }

    /**
     * y-component of velocity
     * 
     * @return y-component of velocity
     */
    public final double getVelocityY() {
        return this.velocityY;
    }

    /**
     * z-component of velocity
     * 
     * @return z-component of velocity
     */
    public final double getVelocityZ() {
        return this.velocityZ;
    }
    // </editor-fold>
    
}
