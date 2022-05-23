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
package de.gnwi.jdpd.samples.harmonicBonds;

import java.util.LinkedList;
import java.util.TreeMap;

/**
 * Harmonic bond chunk where a particle index of all harmonic bonds occurs only 
 * ONCE in the chunk
 * 
 * @author Achim Zielesny
 */
public class HarmonicBondChunk {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Harmonic bond list
     */
    private final LinkedList<HarmonicBond> harmonicBondList;
    
    /**
     * Particle index map
     */
    private final TreeMap<Integer, Integer> particleIndexMap;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public HarmonicBondChunk() {
        this.harmonicBondList = new LinkedList<>();
        this.particleIndexMap = new TreeMap<>();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Returns if bond chunk already contains a particle index
     * 
     * @param aParticleIndex Particle index
     * @return True: Bond chunk already contains particle index, false: Otherwise
     */
    public boolean hasParticleIndex(int aParticleIndex) {
        return this.particleIndexMap.containsKey(aParticleIndex);
    }
    
    /**
     * Adds bond
     * 
     * @param aBond Bond
     * @return True: Operation successful (bond was added to chunk), false: Otherwise
     */
    public boolean addBond(HarmonicBond aBond) {
        if (aBond == null ||
            this.hasParticleIndex(aBond.getIndex1()) ||
            this.hasParticleIndex(aBond.getIndex2())) {
            return false;
        } else {
            this.harmonicBondList.add(aBond);
            this.particleIndexMap.put(aBond.getIndex1(), aBond.getIndex1());
            this.particleIndexMap.put(aBond.getIndex2(), aBond.getIndex2());
            return true;
        }
    }

    /**
     * Returns if harmonic bond chunk has harmonic bonds
     * 
     * @return True: Bond chunk has bonds, false: Otherwise
     */
    public boolean hasBonds() {
        return !this.harmonicBondList.isEmpty();
    }
    
    /**
     * Bond chunk arrays
     * 
     * @return Bond chunk arrays or null if none are available
     */
    public HarmonicBondChunkArrays getBondChunkArrays() {
        if (this.harmonicBondList.isEmpty()) {
            return null;
        } else {
            int[] tmpParticleIndices1 = new int[this.harmonicBondList.size()];
            int[] tmpParticleIndices2 = new int[this.harmonicBondList.size()];
            double[] tmpBondlengths = new double[this.harmonicBondList.size()];
            double[] tmpForceConstants = new double[this.harmonicBondList.size()];
            HarmonicBond.HarmonicBondBehaviour[] tmpHarmonicBondBehaviours = new HarmonicBond.HarmonicBondBehaviour[this.harmonicBondList.size()];
            int tmpIndex = 0;
            for (HarmonicBond tmpBond : this.harmonicBondList) {
                tmpParticleIndices1[tmpIndex] = tmpBond.getIndex1();
                tmpParticleIndices2[tmpIndex] = tmpBond.getIndex2();
                tmpBondlengths[tmpIndex] = tmpBond.getBondLength();
                tmpForceConstants[tmpIndex] = tmpBond.getForceConstant();
                tmpHarmonicBondBehaviours[tmpIndex] = tmpBond.getHarmonicBondBehaviour();
                tmpIndex++;
            }
            return new HarmonicBondChunkArrays(tmpParticleIndices1, tmpParticleIndices2, tmpBondlengths, tmpForceConstants, tmpHarmonicBondBehaviours);
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Returns number of bonds of this bond chunk
     * 
     * @return Number of bonds of this bond chunk
     */
    public int getBondNumber() {
        return this.harmonicBondList.size();
    }
    // </editor-fold>
    
}
