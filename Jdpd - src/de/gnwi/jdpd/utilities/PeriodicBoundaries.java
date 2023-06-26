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
package de.gnwi.jdpd.utilities;

/**
 * Periodic boundaries
 * 
 * @author Achim Zielesny
 */
public class PeriodicBoundaries {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * True: Periodic boundary along x-direction, false: Closed box
     */
    private final boolean isAlongX;

    /**
     * True: Periodic boundary along y-direction, false: Closed box
     */
    private final boolean isAlongY;

    /**
     * True: Periodic boundary along z-direction, false: Closed box
     */
    private final boolean isAlongZ;
    
    /**
     * True: Periodic boundary for at least one direction, false: No periodic 
     * boundaries at all
     */
    private final boolean isPeriodicBoundary;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param anIsAlongX True: Periodic boundary along x-direction, false: Closed box
     * @param anIsAlongY True: Periodic boundary along y-direction, false: Closed box
     * @param anIsAlongZ True: Periodic boundary along z-direction, false: Closed box
     */
    public PeriodicBoundaries(boolean anIsAlongX, boolean anIsAlongY, boolean anIsAlongZ) {
        this.isAlongX = anIsAlongX;
        this.isAlongY = anIsAlongY;
        this.isAlongZ = anIsAlongZ;
        this.isPeriodicBoundary = this.isAlongX || this.isAlongY || this.isAlongZ; 
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * True: Periodic boundary along x-direction, false: Closed box
     * 
     * @return True: Periodic boundary along x-direction, false: Closed box
     */
    public boolean isAlongX() {
        return this.isAlongX;
    }

    /**
     * True: Periodic boundary along y-direction, false: Closed box
     * 
     * @return True: Periodic boundary along y-direction, false: Closed box
     */
    public boolean isAlongY() {
        return this.isAlongY;
    }

    /**
     * True: Periodic boundary along z-direction, false: Closed box
     * 
     * @return True: Periodic boundary along z-direction, false: Closed box
     */
    public boolean isAlongZ() {
        return this.isAlongZ;
    }
    
    /**
     * True: Periodic boundary for at least one direction, false: No periodic 
     * boundaries at all
     * 
     * @return True: Periodic boundary for at least one direction, false: No periodic 
     * boundaries at all
     */
    public boolean isPeriodicBoundary() {
        return this.isPeriodicBoundary;
    }
    // </editor-fold>
    
}
