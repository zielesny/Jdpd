/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2021  Achim Zielesny (achim.zielesny@googlemail.com)
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

/**
 * Harmonic bond chunk related arrays created by HarmonicBondChunk
 * 
 * @author Achim Zielesny
 */
public class HarmonicBondChunkArrays {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Particle indices 1
     */
    private final int[] particleIndices1;

    /**
     * Particle indices 2
     */
    private final int[] particleIndices2;

    /**
     * Bond lengths in DPD units
     */
    private final double[] bondLengths;

    /**
     * Spring force constants in DPD units
     */
    private final double[] forceConstants;
    
    /**
     * Repulsion flags:
     * True: Repulsion for bond is to be calculated, false: Otherwise (no repulsion, attraction only)
     */
    private final HarmonicBond.HarmonicBondBehaviour[] harmonicBondBehaviours;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aParticleIndices1 Particle indices 1
     * @param aParticleIndices2 Particle indices 2
     * @param aBondLengths Bond lengths in DPD units
     * @param aForceConstants Spring force constants in DPD units
     * @param aHarmonicBondBehaviours Behaviours of harmonic bonds
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public HarmonicBondChunkArrays(
        int[] aParticleIndices1,
        int[] aParticleIndices2,
        double[] aBondLengths,
        double[] aForceConstants,
        HarmonicBond.HarmonicBondBehaviour[] aHarmonicBondBehaviours
        ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aParticleIndices1 == null || aParticleIndices1.length == 0) {
            throw new IllegalArgumentException("BondChunkArrays.Constructor: aParticleIndices1 is null/empty.");
        }
        if (aParticleIndices2 == null || aParticleIndices2.length == 0) {
            throw new IllegalArgumentException("BondChunkArrays.Constructor: aParticleIndices2 is null/empty.");
        }
        if (aBondLengths == null || aBondLengths.length == 0) {
            throw new IllegalArgumentException("BondChunkArrays.Constructor: aBondLengths is null/empty.");
        }
        if (aForceConstants == null || aForceConstants.length == 0) {
            throw new IllegalArgumentException("BondChunkArrays.Constructor: aForceConstants is null/empty.");
        }
        if (aParticleIndices1.length != aParticleIndices2.length ||
            aParticleIndices2.length != aBondLengths.length ||
            aBondLengths.length != aForceConstants.length) {
            throw new IllegalArgumentException("BondChunkArrays.Constructor: Arrays have different lengths.");
        }
        // </editor-fold>
        this.particleIndices1 = aParticleIndices1;
        this.particleIndices2 = aParticleIndices2;
        this.bondLengths = aBondLengths;
        this.forceConstants = aForceConstants;
        this.harmonicBondBehaviours = aHarmonicBondBehaviours;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Particle indices 1
     * 
     * @return Particle indices 1
     */
    public int[] getParticleIndices1() {
        return this.particleIndices1;
    }

    /**
     * Particle indices 2
     * 
     * @return Particle indices 2
     */
    public int[] getParticleIndices2() {
        return this.particleIndices2;
    }

    /**
     * Bond lengths in DPD units
     * 
     * @return Bond lengths in DPD units
     */
    public double[] getBondLengths() {
        return this.bondLengths;
    }

    /**
     * Spring force constants in DPD units
     * 
     * @return Spring force constants in DPD units
     */
    public double[] getForceConstants() {
        return this.forceConstants;
    }

    /**
     * Behaviours of harmonic bonds
     * 
     * @return Behaviours of harmonic bonds
     */
    public HarmonicBond.HarmonicBondBehaviour[] getHarmonicBondBehaviours() {
        return this.harmonicBondBehaviours;
    }
    
    /**
     * Returns length of arrays
     * 
     * @return Length of arrays
     */
    public int getLength() {
        return this.particleIndices1.length;
    }
    // </editor-fold>
    
}
