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
package de.gnwi.jdpdsp.particlePosition;

/**
 * Particle position
 * 
 * @author Achim Zielesny
 */
public class ParticlePosition implements Comparable<ParticlePosition> {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * General separator string
     */
    private static final String GENERAL_SEPARATOR = "|";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * x-Position
     */
    private float xPosition;

    /**
     * y-Position
     */
    private float yPosition;

    /**
     * z-Position
     */
    private float zPosition;

    /**
     * Particle token
     */
    private String particleToken;

    /**
     * Moelcule name
     */
    private String moleculeName;

    /**
     * Molecule name plus particle token for comparison
     */
    private String moleculeParticleString;
    
    /**
     * Particle index
     */
    private int particleIndex;
    
    /**
     * Molecule index
     */
    private int moleculeIndex;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     */
    public ParticlePosition() {
        // Do nothing
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Compares moleculeParticleStrings
     *
     * @param anotherParticlePosition ParticlePosition
 instance
     * @return Standard compareTo-result
     */
    @Override
    public int compareTo(ParticlePosition anotherParticlePosition) {
        return this.moleculeParticleString.compareTo(anotherParticlePosition.getMoleculeParticleString());
    }
    
    /**
     * Sets position
     * (No checks are performed)
     * 
     * @param aParticleToken Particle token
     * @param aMoleculeName Molecule name
     * @param aXPosition x-Position
     * @param aYPosition y-Position
     * @param aZPosition z-Position
     * @param aParticleIndex Particle index
     * @param aMoleculeIndex Molecule index
     */
    public final void setPosition(
        String aParticleToken, 
        String aMoleculeName, 
        float aXPosition, 
        float aYPosition, 
        float aZPosition,
        int aParticleIndex,
        int aMoleculeIndex
    ) {
        this.particleToken = aParticleToken;
        this.moleculeName = aMoleculeName;
        this.xPosition = aXPosition;
        this.yPosition = aYPosition;
        this.zPosition = aZPosition;
        this.particleIndex = aParticleIndex;
        this.moleculeIndex = aMoleculeIndex;

        // Compare string is molecule name + particle token with separator
        this.moleculeParticleString = aMoleculeName + GENERAL_SEPARATOR + aParticleToken;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get only)">
    /**
     * Molecule name plus particle token for comparison
     *
     * @return Molecule name plus particle token for comparison
     */
    public String getMoleculeParticleString() {
        return this.moleculeParticleString;
    }

    /**
     * x-Position
     *
     * @return x-Position
     */
    public float getXPosition() {
        return this.xPosition;
    }

    /**
     * y-Position
     *
     * @return y-Position
     */
    public float getYPosition() {
        return this.yPosition;
    }

    /**
     * z-Position
     *
     * @return z-Position
     */
    public float getZPosition() {
        return this.zPosition;
    }

    /**
     * Particle index
     *
     * @return Particle index
     */
    public int getParticleIndex() {
        return this.particleIndex;
    }

    /**
     * Molecule index
     *
     * @return Molecule index
     */
    public int getMoleculeIndex() {
        return this.moleculeIndex;
    }

    /**
     * Particle token
     *
     * @return Particle token
     */
    public String getParticleToken() {
        return this.particleToken;
    }

    /**
     * Molecule name
     *
     * @return Molecule name
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }

	// </editor-fold>
    
}
