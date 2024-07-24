/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2024  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpdsp.utilities;

import org.apache.commons.math3.util.FastMath;

/**
 * Box size
 * 
 * @author Achim Zielesny
 */
public class BoxSize {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Box: x-length
     */
    private final float xLength;
    
    /**
     * Box: y-length
     */
    private final float yLength;

    /**
     * Box: z-length
     */
    private final float zLength;
    
    /**
     * Box: Half x-length
     */
    private final float xHalfLength;
    
    /**
     * Box: Half y-length
     */
    private final float yHalfLength;

    /**
     * Box: Half z-length
     */
    private final float zHalfLength;
    
    /**
     * Box: Double x-length
     */
    private final float xDoubleLength;
    
    /**
     * Box: Double y-length
     */
    private final float yDoubleLength;

    /**
     * Box: Double z-length
     */
    private final float zDoubleLength;
    
    /**
     * Box: Negative half x-length
     */
    private final float negativeXHalfLength;
    
    /**
     * Box: Negative half y-length
     */
    private final float negativeYHalfLength;

    /**
     * Box: Negative half z-length
     */
    private final float negativeZHalfLength;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aXLength Length of box in x-direction
     * @param aYLength Length of box in y-direction
     * @param aZLength Length of box in z-direction
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    public BoxSize (float aXLength, float aYLength, float aZLength) throws IllegalArgumentException {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aXLength <= 0.0f || aYLength <= 0.0f || aZLength <= 0.0f) {
            throw new IllegalArgumentException("BoxSize.Constructor: Length of box is illegal.");
        }
        // </editor-fold>
        this.xLength = aXLength;
        this.xHalfLength = aXLength/2.0f;
        this.xDoubleLength = aXLength + aXLength;
        this.negativeXHalfLength = -this.xHalfLength;
        this.yLength = aYLength;
        this.yHalfLength = aYLength/2.0f;
        this.yDoubleLength = aYLength + aYLength;
        this.negativeYHalfLength = -this.yHalfLength;
        this.zLength = aZLength;
        this.zHalfLength = aZLength/2.0f;
        this.zDoubleLength = aZLength + aZLength;
        this.negativeZHalfLength = -this.zHalfLength;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Box: x-length
     * 
     * @return Box: x-length
     */
    public float getXLength() {
        return this.xLength;
    }
    
    /**
     * Box: y-length
     * 
     * @return Box: y-length
     */
    public float getYLength() {
        return this.yLength;
    }

    /**
     * Box: z-length
     * 
     * @return Box: z-length
     */
    public float getZLength() {
        return this.zLength;
    }
    
    /**
     * Return minimum side length, i.e. minimum of getXLength(), getYLength() 
     * and getZLength()
     * 
     * @return Minimum side length
     */
    public float getMinimumLength() {
        return FastMath.min(this.getXLength(), FastMath.min(this.getYLength(), this.getZLength()));
    }
    
    /**
     * Box: Half x-length
     * 
     * @return Box: Half x-length
     */
    public float getXHalfLength() {
        return this.xHalfLength;
    }
    
    /**
     * Box: Half y-length
     * 
     * @return Box: Half y-length
     */
    public float getYHalfLength() {
        return this.yHalfLength;
    }

    /**
     * Box: Half z-length
     * 
     * @return Box: Half z-length
     */
    public float getZHalfLength() {
        return this.zHalfLength;
    }
    
    /**
     * Box: Double x-length
     * 
     * @return Box: Double x-length
     */
    public float getXDoubleLength() {
        return this.xDoubleLength;
    }
    
    /**
     * Box: Double y-length
     * 
     * @return Box: Double y-length
     */
    public float getYDoubleLength() {
        return this.yDoubleLength;
    }

    /**
     * Box: Double z-length
     * 
     * @return Box: Double z-length
     */
    public float getZDoubleLength() {
        return this.zDoubleLength;
    }
    
    /**
     * Box: Negative half x-length
     * 
     * @return Box: Negative half x-length
     */
    public float getNegativeXHalfLength() {
        return this.negativeXHalfLength;
    }
    
    /**
     * Box: Negative half y-length
     * 
     * @return Box: Negative half y-length
     */
    public float getNegativeYHalfLength() {
        return this.negativeYHalfLength;
    }

    /**
     * Box: Negative half z-length
     * 
     * @return Box: Negative half z-length
     */
    public float getNegativeZHalfLength() {
        return this.negativeZHalfLength;
    }
    // </editor-fold>
    
}
