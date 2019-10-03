/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2019  Achim Zielesny (achim.zielesny@googlemail.com)
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
 * Molecule fixation description
 * 
 * @author Achim Zielesny
 */
public class MoleculeFixationDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;

    /**
     * True: x-component of all particles of molecule is fixed, false: Otherwise
     */
    private final boolean isFixedX;

    /**
     * True: y-component of all particles of molecule is fixed, false: Otherwise
     */
    private final boolean isFixedY;
    
    /**
     * True: z-component of all particles of molecule is fixed, false: Otherwise
     */
    private final boolean isFixedZ;

    /**
     * Maximum time step for application
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
     * @param anIsFixedX True: x-component of all particles of molecule is fixed, false: Otherwise
     * @param anIsFixedY True: y-component of all particles of molecule is fixed, false: Otherwise
     * @param anIsFixedZ True: z-component of all particles of molecule is fixed, false: Otherwise
     * @param aMaxTimeStep Maximum time step for application
     */
    public MoleculeFixationDescription(
        String aMoleculeName,
        boolean anIsFixedX,
        boolean anIsFixedY,
        boolean anIsFixedZ,
        int aMaxTimeStep
        ) {
        this.moleculeName = aMoleculeName;
        this.isFixedX = anIsFixedX;
        this.isFixedY = anIsFixedY;
        this.isFixedZ = anIsFixedZ;
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
     * True: x-component of all particles of molecule is fixed, false: Otherwise
     * 
     * @return True: x-component of all particles of molecule is fixed, false: Otherwise
     */
    public boolean isFixedX() {
        return this.isFixedX;
    }

    /**
     * True: y-component of all particles of molecule is fixed, false: Otherwise
     * 
     * @return True: y-component of all particles of molecule is fixed, false: Otherwise
     */
    public boolean isFixedY() {
        return this.isFixedY;
    }

    /**
     * True: z-component of all particles of molecule is fixed, false: Otherwise
     * 
     * @return True: z-component of all particles of molecule is fixed, false: Otherwise
     */
    public boolean isFixedZ() {
        return this.isFixedZ;
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
