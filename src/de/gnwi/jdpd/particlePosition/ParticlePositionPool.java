/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2022  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.particlePosition;

import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Pool for particle positions
 * 
 * @author Achim Zielesny
 */
public class ParticlePositionPool {
    
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Particle position queue
     */
    private final ConcurrentLinkedQueue<ParticlePosition> particlePositionQueue;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public ParticlePositionPool() {
        this.particlePositionQueue = new ConcurrentLinkedQueue<>();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Returns particle position
     * 
     * @return Particle position
     */
    public ParticlePosition getParticlePosition() {
        ParticlePosition tmpParticlePosition = this.particlePositionQueue.poll();
        if (tmpParticlePosition == null) {
            return new ParticlePosition();
        } else {
            return tmpParticlePosition;
        }
    }
    
    /**
     * Set particle position for reuse
     * 
     * @param aParticlePosition Particle position
     */
    public void setParticlePositionForReuse(ParticlePosition aParticlePosition) {
        if (aParticlePosition != null) {
            this.particlePositionQueue.add(aParticlePosition);
        }
    }
    
    /**
     * Clear pool
     */
    public void clear() {
        this.particlePositionQueue.clear();
    }
    // </editor-fold>
    
}
