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
package de.gnwi.jdpd.parameters;

import de.gnwi.jdpd.utilities.Electrostatics;
import de.gnwi.jdpd.utilities.GravitationalAcceleration;
import org.apache.commons.math3.util.FastMath;

/**
 * Simulation interaction parameters
 * 
 * @author Achim Zielesny
 */
public class InteractionDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     */
    private final double temperature;

    /**
     * DPD sigma parameter in DPD units
     */
    private final double dpdSigma;

    /**
     * DPD sigma parameter divided by root of time step length in DPD units
     */
    private final double dpdSigmaDivRootTimeStepLength;

    /**
     * DPD gamma parameter in DPD units
     */
    private final double dpdGamma;
    
    /**
     * Flag for use of Gaussian random for random force
     * 
     * True : Random DPD force is driven by random variable with 
     *        Gaussian distribution with zero mean and unit variance
     *        (slower)
     * False: Random DPD force is driven by random variable with uniform
     *        distribution with zero mean and unit variance (faster)
     */    
    private final boolean isGaussianRandomDpdForce;
    
    /**
     * Number of particle types
     */
    private final int particleTypeNumber;
    
    /**
     * Conservative force parameters a(ij)
     * a(ij) = aij[i * particleTypeNumber + j]
     */
    private final double[][] aij;
    
    /**
     * Electrostatics parameters
     */
    private final Electrostatics electrostatics;
    
    /**
     * Gravitational acceleration
     */
    private final GravitationalAcceleration gravitationalAcceleration;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables for PNHLN integration type">
    /**
     * Ksi value for PNHLN integration type
     */
    private double pnhln_Ksi;
    
    /**
     * True: G value for PNHLN integration type is calculated, false: Otherwise
     */
    private boolean isPnhln_G_Calculation;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aTemperature Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * @param aDpdSigma DPD sigma parameter in DPD units
     * @param anIsGaussianRandomDpdForce Flag for use of Gaussian random for random force
     * @param aParticleTypeNumber Number of particle types
     * @param anAij Conservative force repulsion parameters a(ij)
     * @param anElectrostatics Electrostatics parameters (may be null)
     * @param aGravitationalAcceleration Gravitational acceleration
     * @param aSimulationDescription Simulation description
     */
    public InteractionDescription(
        double aTemperature,
        double aDpdSigma,
        boolean anIsGaussianRandomDpdForce,
        int aParticleTypeNumber,
        double[][] anAij,
        Electrostatics anElectrostatics,
        GravitationalAcceleration aGravitationalAcceleration,
        SimulationDescription aSimulationDescription
    ) {
        this.temperature = aTemperature;
        this.dpdSigma = aDpdSigma;
        // DPD sigma parameter divided by root of time step length in DPD units
        this.dpdSigmaDivRootTimeStepLength = this.dpdSigma/FastMath.sqrt(aSimulationDescription.getTimeStepLength());
        this.isGaussianRandomDpdForce = anIsGaussianRandomDpdForce;
        this.particleTypeNumber = aParticleTypeNumber;
        this.aij = anAij;
        this.electrostatics = anElectrostatics;
        this.gravitationalAcceleration = aGravitationalAcceleration;
        // DPD gamma parameter
        // NOTE: this.parameters.getInteractionDescription().temperature is defined in kT fraction (k: Boltzmann constant)
        this.dpdGamma = this.dpdSigma * this.dpdSigma / (2.0 * this.temperature);
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties for PNHLN integration type (get/set)">
    /**
     * Ksi value for PNHLN integration type
     * 
     * @return Ksi value for PNHLN integration type
     */
    public double getPnhln_Ksi() {
        return this.pnhln_Ksi;
    }

    /**
     * Ksi value for PNHLN integration type
     * 
     * @param aValue Ksi value for PNHLN integration type
     */
    public void setPnhln_Ksi(double aValue) {
        this.pnhln_Ksi = aValue;
    }
    
    /**
     * True: G value for PNHLN integration type is calculated, false: Otherwise
     * 
     * @return True: G value for PNHLN integration type is calculated, false: Otherwise
     */
    public boolean isPnhln_G_Calculation() {
        return this.isPnhln_G_Calculation;
    }
    
    /**
     * True: G value for PNHLN integration type is calculated, false: Otherwise
     * 
     * @param aValue G value calculation flag for PNHLN integration type
     */
    public void setPnhln_G_Calculation(boolean aValue) {
        this.isPnhln_G_Calculation = aValue;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * 
     * @return Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     */
    public double getTemperature() {
        return this.temperature;
    }

    /**
     * DPD sigma parameter in DPD units
     * 
     * @return DPD sigma parameter in DPD units
     */
    public double getDpdSigma() {
        return this.dpdSigma;
    }

    /**
     * DPD sigma parameter divided by root of time step length in DPD units
     * 
     * @return DPD sigma parameter divided by root of time step length in DPD units
     */
    public double getDpdSigmaDivRootTimeStepLength() {
        return this.dpdSigmaDivRootTimeStepLength;
    }

    /**
     * Flag for use of Gaussian random for random force
     * 
     * @return True : Random DPD force is driven by random variable with 
     *                Gaussian distribution with zero mean and unit variance
     *                (slower)
     *         False: Random DPD force is driven by random variable with uniform
     *                distribution with zero mean and unit variance (faster)
     */
    public boolean isGaussianRandomDpdForce() {
        return this.isGaussianRandomDpdForce;
    }
    
    /**
     * Number of particle types
     * 
     * @return Number of particle types
     */
    public int getParticleTypeNumber() {
        return this.particleTypeNumber;
    }
    
    /**
     * Conservative force parameters a(ij)
     * 
     * @return Conservative force parameters a(ij)
     */
    public double[][] getAij() {
        return this.aij;
    }

    /**
     * Electrostatics parameters
     * 
     * @return Electrostatics parameters or null if none are available
     */
    public Electrostatics getElectrostatics() {
        return this.electrostatics;
    }
    
    /**
     * Returns if electrostatics parameters are available
     * 
     * @return True: Electrostatics parameters are available, false: Otherwise
     */
    public boolean hasElectrostatics() {
        return this.electrostatics != null;
    }

    /**
     * Gravitational acceleration
     * 
     * @return Gravitational acceleration
     */
    public GravitationalAcceleration getGravitationalAcceleration() {
        return this.gravitationalAcceleration;
    }

    /**
     * True: Gravitational acceleration is defined, false: Otherwise
     * 
     * @return True: Gravitational acceleration is defined, false: Otherwise
     */
    public boolean isGravitationalAcceleration() {
        return this.gravitationalAcceleration.isGravitationalAcceleration();
    }
    
    /**
     * DPD gamma parameter in DPD units
     * 
     * @return DPD gamma parameter in DPD units
     */
    public double getDpdGamma() {
        return this.dpdGamma;
    }
    // </editor-fold>
    
}
