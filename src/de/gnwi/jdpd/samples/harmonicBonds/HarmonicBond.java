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

import java.util.Arrays;

/**
 * Harmonic bond
 * 
 * @author Achim Zielesny
 */
public class HarmonicBond {

    // <editor-fold defaultstate="collapsed" desc="Public enums">
    /**
     * Behaviour of harmonic bond
     */
    public enum HarmonicBondBehaviour {
        
        /**
         * Harmonic bond is repulsive and attractive (default)
         */
        DEFAULT,
        /**
         * Harmonic bond is only repulsive if squeezed (but NOT attractive if
         * elongated)
         */
        REPULSIVE,
        /**
         * Harmonic bond is only attractive if elongated (but NOT repulsive if
         * squeezed)
         */
        ATTRACTIVE;

        /**
         * Harmonic bond behaviour representations
         * 
         * @return Harmonic bond behaviour representations
         */
        public static String[] getHarmonicBondBehaviourRepresentations() {
            String[] tmpHarmonicBondBehaviourRepresentations = new String[HarmonicBondBehaviour.values().length];
            int tmpIndex = 0;
            for (HarmonicBondBehaviour tmpHarmonicBondBehaviour : HarmonicBondBehaviour.values()) {
                tmpHarmonicBondBehaviourRepresentations[tmpIndex++] = tmpHarmonicBondBehaviour.toString();
            }
            Arrays.sort(tmpHarmonicBondBehaviourRepresentations);
            return tmpHarmonicBondBehaviourRepresentations;
        }

        /**
         * Default harmonic bond behaviour representation
         * 
         * @return Default harmonic bond behaviour representation
         */
        public static String getDefaultHarmonicBondBehaviourRepresentation() {
            return HarmonicBondBehaviour.getDefaultHarmonicBondBehaviour().toString();
        }

        /**
         * Default harmonic bond behaviour 
         * 
         * @return Default harmonic bond behaviour 
         */
        public static HarmonicBondBehaviour getDefaultHarmonicBondBehaviour() {
            return HarmonicBondBehaviour.DEFAULT;
        }
        
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Index 1
     */
    private final int index1;

    /**
     * Index 2
     */
    private final int index2;
    
    /**
     * Bond length in DPD units
     */
    private final double bondLength;

    /**
     * Spring force constant in DPD units
     */
    private final double forceConstant;
    
    /**
     * Behaviour of harmonic bond
     */
    private final HarmonicBondBehaviour harmonicBondBehaviour;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param anIndex1 Index of particle 1
     * @param anIndex2 Index of particle 2
     * @param aBondLength Bond length in DPD units
     * @param aForceConstant Spring force constant in DPD units
     * @param aHarmonicBondBehaviour Behaviour of harmonic bond
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public HarmonicBond(
        int anIndex1,
        int anIndex2,
        double aBondLength,
        double aForceConstant,
        HarmonicBondBehaviour aHarmonicBondBehaviour
        ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anIndex1 < 0) {
            throw new IllegalArgumentException("HarmonicBond.Constructor: anIndex1 < 0.");
        }
        if (anIndex2 < 0) {
            throw new IllegalArgumentException("HarmonicBond.Constructor: anIndex2 < 0.");
        }
        if (aBondLength <= 0.0) {
            throw new IllegalArgumentException("HarmonicBond.Constructor: aBondLength <= 0.0.");
        }
        if (aForceConstant <= 0.0) {
            throw new IllegalArgumentException("HarmonicBond.Constructor: aForceConstant <= 0.0.");
        }
        // </editor-fold>
        this.index1 = anIndex1;
        this.index2 = anIndex2;
        this.bondLength = aBondLength;
        this.forceConstant = aForceConstant;
        this.harmonicBondBehaviour = aHarmonicBondBehaviour;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Index 1
     * 
     * @return Index 1
     */
    public int getIndex1() {
        return this.index1;
    }

    /**
     * Index 2
     * 
     * @return Index 2
     */
    public int getIndex2() {
        return this.index2;
    }

    /**
     * Bond length in DPD units
     * 
     * @return Bond length in DPD units
     */
    public double getBondLength() {
        return this.bondLength;
    }

    /**
     * Spring force constant in DPD units
     * 
     * @return Spring force constant in DPD units
     */
    public double getForceConstant() {
        return this.forceConstant;
    }
    
    /**
     * Behaviour of harmonic bond
     * 
     * @return Behaviour of harmonic bond
     */
    public HarmonicBondBehaviour getHarmonicBondBehaviour() {
        return this.harmonicBondBehaviour;
    }
    // </editor-fold>
    
}
