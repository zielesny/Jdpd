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
package de.gnwi.jdpd.samples.harmonicBonds;

/**
 * Particle-pair harmonic bond
 * 
 * @author Achim Zielesny
 */
public class ParticlePairHarmonicBond {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Particle-pair key
     */
    private final String particlePairKey;

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
    private final HarmonicBond.HarmonicBondBehaviour harmonicBondBehaviour;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aParticlePairKey Particle-pair key
     * @param aBondLength Bond length in DPD units
     * @param aForceConstant Spring force constant in DPD units
     * @param aHarmonicBondBehaviour Behaviour of harmonic bond
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairHarmonicBond(
        String aParticlePairKey,
        double aBondLength,
        double aForceConstant,
        HarmonicBond.HarmonicBondBehaviour aHarmonicBondBehaviour
        ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aParticlePairKey == null || aParticlePairKey.isEmpty()) {
            throw new IllegalArgumentException("Bond.Constructor: aParticlePairKey is null/empty.");
        }
        if (aBondLength <= 0.0) {
            throw new IllegalArgumentException("Bond.Constructor: aBondLength <= 0.0.");
        }
        if (aForceConstant <= 0.0) {
            throw new IllegalArgumentException("Bond.Constructor: aForceConstant <= 0.0.");
        }
        // </editor-fold>
        this.particlePairKey = aParticlePairKey;
        this.bondLength = aBondLength;
        this.forceConstant = aForceConstant;
        this.harmonicBondBehaviour = aHarmonicBondBehaviour;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Particle-pair key
     * 
     * @return Particle-pair key
     */
    public String getParticlePairKey() {
        return this.particlePairKey;
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
    public HarmonicBond.HarmonicBondBehaviour getHarmonicBondBehaviour() {
        return this.harmonicBondBehaviour;
    }
    // </editor-fold>
    
}
