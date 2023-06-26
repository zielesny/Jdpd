/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpdsp.parameters;

import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.DoubleAdder;

/**
 * Objects for testing purposes
 * 
 * @author Achim Zielesny
 */
public class TestObjects {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Queue of particle pair arrays [length = 2]
     */
    private final ConcurrentLinkedQueue<int[]> particleIndexPairQueue;

    /**
     * Queue of particle pair random strings
     */
    private final ConcurrentLinkedQueue<String> particleIndexPairRandomQueue;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Particle index pair counter
     */
    private AtomicInteger particleIndexPairCounter;

    /**
     * Double adder
     */
    private DoubleAdder doubleAdder;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aParticleIndexPairQueue Queue of particle pair arrays [length = 2]
     * @param aParticleIndexPairRandomQueue Queue of particle pair random strings
     * @param aParticleIndexPairNumber Particle index pair counter
     * @param aDoubleAdder Double adder
     */
    public TestObjects(
        ConcurrentLinkedQueue<int[]> aParticleIndexPairQueue,
        ConcurrentLinkedQueue<String> aParticleIndexPairRandomQueue,
        AtomicInteger aParticleIndexPairNumber,
        DoubleAdder aDoubleAdder
    ) {
        this.particleIndexPairQueue = aParticleIndexPairQueue;
        this.particleIndexPairRandomQueue = aParticleIndexPairRandomQueue;
        this.particleIndexPairCounter = aParticleIndexPairNumber;
        this.doubleAdder = aDoubleAdder;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Queue of particle pair arrays [length = 2]
     * 
     * @return Queue of particle pair arrays [length = 2]
     */
    public ConcurrentLinkedQueue<int[]> getParticleIndexPairQueue() {
        return this.particleIndexPairQueue;
    }

    /**
     * Queue of particle pair random strings
     * 
     * @return Queue of particle pair random strings
     */
    public ConcurrentLinkedQueue<String> getParticleIndexPairRandomQueue() {
        return this.particleIndexPairRandomQueue;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get/set)">
    /**
     * Particle index pair counter
     * 
     * @return Particle index pair counter
     */
    public AtomicInteger getParticleIndexPairCounter() {
        return this.particleIndexPairCounter;
    }
    
    /**
     * Particle index pair counter
     * 
     * @param aParticleIndexPairCounter Particle index pair counter
     */
    public void setParticleIndexPairCounter(AtomicInteger aParticleIndexPairCounter) {
        this.particleIndexPairCounter = aParticleIndexPairCounter;
    }

    /**
     * Double adder
     * 
     * @return Double adder
     */
    public DoubleAdder getDoubleAdder() {
        return this.doubleAdder;
    }
    
    /**
     * Double adder
     * 
     * @param aDoubleAdder Double adder
     */
    public void setDoubleAdder(DoubleAdder aDoubleAdder) {
        this.doubleAdder = aDoubleAdder;
    }
    // </editor-fold>
    
}
