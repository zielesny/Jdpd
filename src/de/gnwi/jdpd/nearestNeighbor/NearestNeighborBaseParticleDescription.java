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
package de.gnwi.jdpd.nearestNeighbor;

import java.util.LinkedList;

/**
 * Nearest-neighbor base particle description (i.e. particles for which a 
 * non-bonded nearest-neighbor particle of another molecule is to be determined)
 * 
 * @author Achim Zielesny
 */
public class NearestNeighborBaseParticleDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule of base particle
     */
    private final String moleculeName;
    
    /**
     * Particle token of base particle
     */
    private final String baseParticleToken;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * List of base particle indices in particle arrays
     */
    private LinkedList<Integer> baseParticleIndexList;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aMoleculeName Name of molecule
     * @param aBaseParticleToken Base particle token
     */
    public NearestNeighborBaseParticleDescription(
        String aBaseParticleToken,
        String aMoleculeName
    ) {
        this.moleculeName = aMoleculeName;
        this.baseParticleToken = aBaseParticleToken;
        this.baseParticleIndexList = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Name of molecule of base particle
     * 
     * @return Name of molecule of base particle
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }

    /**
     * Particle token of base particle
     * 
     * @return Particle token of base particle
     */
    public String getBaseParticleToken() {
        return this.baseParticleToken;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get/set)">
    /**
     * List of base particle indices in particle arrays
     * 
     * @return List of base particle indices in particle arrays
     */
    public LinkedList<Integer> getBaseParticleIndexList() {
        return this.baseParticleIndexList;
    }

    /**
     * List of base particle indices in particle arrays
     * 
     * @param aBaseParticleIndexList List of base particle indices in particle arrays
     */
    public void setBaseParticleIndexList(LinkedList<Integer> aBaseParticleIndexList) {
        this.baseParticleIndexList = aBaseParticleIndexList;
    }
    // </editor-fold>
    
}
