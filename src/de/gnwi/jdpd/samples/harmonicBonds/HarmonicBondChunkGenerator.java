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

/**
 * Generator for parallelisation safe harmonic bond chunks
 * 
 * @author Achim Zielesny
 */
public class HarmonicBondChunkGenerator {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    private static final String PLUS = " + ";
    private static final String EQUALS = " = ";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Bond chunk list
     */
    private final LinkedList<HarmonicBondChunk> bondChunkList;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public HarmonicBondChunkGenerator() {
        this.bondChunkList = new LinkedList<>();
        this.bondChunkList.add(new HarmonicBondChunk());
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Adds bond (may be a doublette)
     * 
     * @param aBond Bond
     * @return True: Operation successful (bond was added), false: Otherwise
     */
    public boolean addBond(HarmonicBond aBond) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aBond == null) {
            return false;
        }
        // </editor-fold>
        boolean tmpIsAddedToBondChunk = false;
        for (HarmonicBondChunk tmpBondChunk : this.bondChunkList) {
            if (tmpBondChunk.addBond(aBond)) {
                tmpIsAddedToBondChunk = true;
                break;
            }
        }
        if (!tmpIsAddedToBondChunk) {
            HarmonicBondChunk tmpBondChunk = new HarmonicBondChunk();
            tmpBondChunk.addBond(aBond);
            this.bondChunkList.add(tmpBondChunk);
        }
        return true;
    }
    
    /**
     * Returns if a bond chunk has bonds
     * 
     * @return True: Bond chunk has bonds, false: Otherwise
     */
    public boolean hasBonds() {
        for (HarmonicBondChunk tmpBondChunk : this.bondChunkList) {
            if (tmpBondChunk.hasBonds()) {
                return true;
            }
        }
        return false;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Returns list of bond chunk arrays
     * @return List of bond chunk arrays or null if none are available
     */
    public LinkedList<HarmonicBondChunkArrays> getBondChunkArraysList() {
        if (this.hasBonds()) {
            LinkedList<HarmonicBondChunkArrays> tmpBondChunkArraysList = new LinkedList<>();
            for (HarmonicBondChunk tmpBondChunk : this.bondChunkList) {
                tmpBondChunkArraysList.add(tmpBondChunk.getBondChunkArrays());
            }
            return tmpBondChunkArraysList;
        } else {
            return null;
        }
    }
    
    /**
     * Returns number of bond chunks
     * 
     * @return Number of bond chunks (0 if no bonds were added)
     */
    public int getBondChunkNumber() {
        if (this.hasBonds()) {
            return this.bondChunkList.size();
        } else {
            return 0;
        }
    }
    
    /**
     * Returns number of bonds of all bond chunks
     * 
     * @return Number of bonds of all bond chunks
     */
    public int getBondNumber() {
        if (this.hasBonds()) {
            int tmpBondNumber = 0;
            for (HarmonicBondChunk tmpBondChunk : this.bondChunkList) {
                tmpBondNumber += tmpBondChunk.getBondNumber();
            }
            return tmpBondNumber;
        } else {
            return 0;
        }
    }
    
    /**
     * Returns bond number info
     * 
     * @return Bond number info or empty string if no bond exists.
     */
    public String getBondInfo() {
        if (this.hasBonds()) {
            // Arbitrary but sufficient capacity in most cases
            StringBuilder tmpBuffer = new StringBuilder(1000);
            int tmpBondNumber = 0;
            for (HarmonicBondChunk tmpBondChunk : this.bondChunkList) {
                if (tmpBuffer.length() > 0) {
                    tmpBuffer.append(PLUS);
                }
                tmpBuffer.append(String.valueOf(tmpBondChunk.getBondNumber()));
                tmpBondNumber += tmpBondChunk.getBondNumber();
            }
            tmpBuffer.append(EQUALS);
            tmpBuffer.append(String.valueOf(tmpBondNumber));
            return tmpBuffer.toString();
        } else {
            return "";
        }
    }
    // </editor-fold>
    
}
