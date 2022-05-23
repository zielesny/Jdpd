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
package de.gnwi.jdpd.utilities;

/**
 * Electrostatics parameters
 * 
 * @author Achim Zielesny
 */
public class Electrostatics {

    // <editor-fold defaultstate="collapsed" desc="Public enums">
    /**
     * Type of charge distribution
     */
    public enum ChargeDistributionType {
        
        /**
         * No charge distribution, i.e. point charges
         */
        NONE,
        /**
         * Alejandre charge distribution
         * 
         * Source: 
         * M. González-Melchor, E. Mayoral, M. E. Velázquez and J. Alejandre
         * Electrostatic interactions in dissipative particle dynamics using 
         * the Ewald sums
         * The Journal of Chemical Physics 125, 2006, 224107
         */
        ALEJANDRE
        
    }

    /**
     * Splitting type
     */
    public enum SplittingType {
        
        /**
         * No splitting
         */
        NONE,
        /**
         * Fanourgakis splitting
         * 
         * Source: 
         * G. S. Fanourgakis
         * An Extension of Wolf’s Method for the Treatment of Electrostatic 
         * Interactions: Application to Liquid Water and Aqueous Solutions
         * J. Phys. Chem. B 119, 2015, 1974−1985
         */
        FANOURGAKIS
        
    }
    // </editor-fold>
    //
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
    
    /**
     * Electrostatics coupling constant
     */
    private final double electrostaticsCoupling;
    
    /**
     * Charge distribution type
     */
    private final ChargeDistributionType chargeDistributionType;
    
    /**
     * Alejandre decay length
     */
    private final double decayLengthAlejandre;
    
    /**
     * Splitting type
     */
    private final SplittingType splittingType;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * (No checks are performed)
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

        this.electrostaticsCoupling = 1;
        this.chargeDistributionType = ChargeDistributionType.NONE;
        this.decayLengthAlejandre = 1;
        this.splittingType = SplittingType.NONE;
    }

    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aCutOffLength Cut-off length in DPD units
     * @param aMaximumAbsoluteForceValue Maximum absolute force value in DPD units
     * @param anEffectiveExponent Effective exponent
     * @param aDampingDistance Damping distance in DPD units
     * @param aDampingFactor Damping factor
     * @param anElectrostaticsCoupling Electrostatics coupling constant
     * @param aChargeDistributionType Charge distribution type
     * @param aDecayLengthAlejandre Alejandre decay length
     * @param aSplittingType Splitting type
     */
    public Electrostatics(
        double aCutOffLength,
        double aMaximumAbsoluteForceValue,
        double anEffectiveExponent,
        double aDampingDistance,
        double aDampingFactor,
        double anElectrostaticsCoupling,
        ChargeDistributionType aChargeDistributionType,
        double aDecayLengthAlejandre,
        SplittingType aSplittingType
        ) {
        this.cutOffLength = aCutOffLength;
        this.maximumAbsoluteForceValue = aMaximumAbsoluteForceValue;
        this.effectiveExponent = anEffectiveExponent;
        this.dampingDistance = aDampingDistance;
        this.dampingFactor = aDampingFactor;
        this.electrostaticsCoupling = anElectrostaticsCoupling;
        this.chargeDistributionType = aChargeDistributionType;
        this.decayLengthAlejandre = aDecayLengthAlejandre;
        this.splittingType = aSplittingType;

        this.effectiveChargeFactor = 1;
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
   
    /**
     * Electrostatics coupling constant
     * 
     * @return Electrostatics coupling constant
     */
    public double getElectrostaticsCoupling() {
        return this.electrostaticsCoupling;
    }
    
    /**
     * Charge distribution type
     * 
     * @return Charge distribution type
     */
    public ChargeDistributionType getChargeDistributionType() {
        return this.chargeDistributionType;
    }
    
    /**
     * Alejandre decay length
     * 
     * @return Alejandre decay length
     */
    public double getDecayLengthAlejandre() {
        return this.decayLengthAlejandre;
    }
    
    /**
     * Splitting type
     * 
     * @return Splitting type
     */
    public SplittingType getSplittingType() {
        return this.splittingType;
    }
    // </editor-fold>
    
}
