/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpd.parameters;

/**
 * Simulation counts parameters
 * 
 * @author Achim Zielesny
 */
public class SimulationCounts {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Number of particles in simulation box
     */
    private final int particleNumber;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aParticleNumber Number of particles in simulation box
     */
    public SimulationCounts(
        int aParticleNumber) {
        this.particleNumber = aParticleNumber;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Number of particles in simulation box
     * 
     * @return Number of particles in simulation box
     */
    public int getParticleNumber() {
        return this.particleNumber;
    }
    // </editor-fold>
        
}
