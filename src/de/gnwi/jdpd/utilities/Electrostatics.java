/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2019  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.utilities;

/**
 * Electrostatics parameters
 * 
 * @author Achim Zielesny
 */
public class Electrostatics {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Cut-off length in DPD units
     */
    private final double cutOffLength;

    /**
     * Maximum absolute force value in DPD units
     */
    private final double maximumAbsoluteForceValue;

    /**
     * Effective charge factor
     */
    private final double effectiveChargeFactor;

    /**
     * Effective exponent
     */
    private final double effectiveExponent;

    /**
     * Damping distance in DPD units
     */
    private final double dampingDistance;

    /**
     * Damping factor
     */
    private final double dampingFactor;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * NOTE: NO checks are performed.
     * 
     * @param aCutOffLength Cut-off length in DPD units
     * @param aMaximumAbsoluteForceValue Maximum absolute force value in DPD units
     * @param anEffectiveChargeFactor Effective charge factor
     * @param anEffectiveExponent Effective exponent
     * @param aDampingDistance Damping distance in DPD units
     * @param aDampingFactor Damping factor
     */
    public Electrostatics(
        double aCutOffLength,
        double aMaximumAbsoluteForceValue,
        double anEffectiveChargeFactor,
        double anEffectiveExponent,
        double aDampingDistance,
        double aDampingFactor
        ) {
        this.cutOffLength = aCutOffLength;
        this.maximumAbsoluteForceValue = aMaximumAbsoluteForceValue;
        this.effectiveChargeFactor = anEffectiveChargeFactor;
        this.effectiveExponent = anEffectiveExponent;
        this.dampingDistance = aDampingDistance;
        this.dampingFactor = aDampingFactor;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Cut-off length in DPD units
     * 
     * @return Cut-off length in DPD units
     */
    public double getCutOffLength() {
        return this.cutOffLength;
    }

    /**
     * Maximum absolute force value in DPD units
     * 
     * @return Maximum absolute force value in DPD units
     */
    public double getMaximumAbsoluteForceValue() {
        return this.maximumAbsoluteForceValue;
    }

    /**
     * Effective charge factor
     * 
     * @return Effective charge factor
     */
    public double getEffectiveChargeFactor() {
        return this.effectiveChargeFactor;
    }

    /**
     * Effective exponent
     * 
     * @return Effective exponent
     */
    public double getEffectiveExponent() {
        return this.effectiveExponent;
    }

    /**
     * Damping distance in DPD units
     * 
     * @return Damping distance in DPD units
     */
    public double getDampingDistance() {
        return this.dampingDistance;
    }

    /**
     * Damping factor
     * 
     * @return Damping factor
     */
    public double getDampingFactor() {
        return this.dampingFactor;
    }
    // </editor-fold>
    
}
