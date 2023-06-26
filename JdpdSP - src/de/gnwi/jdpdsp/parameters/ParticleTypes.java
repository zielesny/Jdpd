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

import java.util.HashMap;

/**
 * Particle types
 *
 * @author Achim Zielesny
 */
public class ParticleTypes {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Number of particle types
     */
    private final int particleTypeNumber;
    
    /**
     * Map of particle token to particle index. This index refers to particle related arrays.
     */
    private final HashMap<String, Integer> particleTokenToIndexMap;
    
    /**
     * Tokens of particles (e.g. "H2O")
     */
    private final String[] particleTokens;
    
    /**
     * Charges of particles
     */
    private final float[] charges;
    
    /**
     * Molar masses of particles in g/mol
     */
    private final float[] molarMasses;
    
    /**
     * Minimum molar mass in g/mol
     */
    private final float minimumMolarMass;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)

     * @param aParticleTokenToIndexMap Map of particle token to particle index. This index refers to particle related arrays.
     * @param aParticleTokens Tokens of particles (e.g. "H2O")
     * @param aCharges Charges of particles
     * @param aMolarMasses Molar masses of particles in g/mol
     */
    public ParticleTypes(
        HashMap<String, Integer> aParticleTokenToIndexMap,
        String[] aParticleTokens,
        float[] aCharges,
        float[] aMolarMasses) {
        this.particleTypeNumber = aParticleTokenToIndexMap.size();
        this.particleTokenToIndexMap = aParticleTokenToIndexMap;
        this.particleTokens = aParticleTokens;
        this.charges = aCharges;
        this.molarMasses = aMolarMasses;
        // Set minimum molar mass
        float tmpMinimumMolarMass = aMolarMasses[0];
        for (int i = 1; i < aMolarMasses.length; i++) {
            if (aMolarMasses[i] < tmpMinimumMolarMass) {
                tmpMinimumMolarMass = aMolarMasses[i];
            }
        }
        this.minimumMolarMass = tmpMinimumMolarMass;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Number of particle types
     * 
     * @return Number of particle types
     */
    public int getParticleTypeNumber() {
        return this.particleTypeNumber;
    }
    
    /**
     * Index of particle token
     * 
     * @param aToken Particle token
     * @return Index of particle token
     */
    public int getIndex(String aToken) {
        return this.particleTokenToIndexMap.get(aToken);
    }
    
    /**
     * Tokens of particle (e.g. "H2O")
     * 
     * @param anIndex Particle index
     * @return Tokens of particle (e.g. "H2O")
     */
    public String getParticleToken(int anIndex) {
        return this.particleTokens[anIndex];
    }
    
    /**
     * Molar mass of particle in g/mol
     * 
     * @param aToken Particle token
     * @return Molar mass of particle in g/mol
     */
    public float getMolarMass(String aToken) {
        return this.molarMasses[this.particleTokenToIndexMap.get(aToken)];
    }
    
    /**
     * Molar mass of particle in g/mol
     * 
     * @param anIndex Particle index
     * @return Molar mass of particle in g/mol
     */
    public float getMolarMass(int anIndex) {
        return this.molarMasses[anIndex];
    }
    
    /**
     * Minimum molar mass in g/mol
     * 
     * @return Minimum molar mass in g/mol
     */
    public float getMinimumMolarMass() {
        return this.minimumMolarMass;
    }
    
    /**
     * Charge of particle
     * 
     * @param aToken Particle token
     * @return Charge of particle
     */
    public float getCharge(String aToken) {
        return this.charges[this.particleTokenToIndexMap.get(aToken)];
    }
    
    /**
     * Charge of particle
     * 
     * @param anIndex Particle index
     * @return Charge of particle
     */
    public float getCharge(int anIndex) {
        return this.charges[anIndex];
    }
    // </editor-fold>
    
}
