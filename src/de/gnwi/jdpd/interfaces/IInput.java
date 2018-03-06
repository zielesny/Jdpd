/**
 * Jdpd - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
 * Copyright (C) 2018  Achim Zielesny (achim.zielesny@googlemail.com)
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
package de.gnwi.jdpd.interfaces;

import de.gnwi.jdpd.movement.MoleculeAccelerationDescription;
import de.gnwi.jdpd.movement.MoleculeBoundaryDescription;
import de.gnwi.jdpd.utilities.Factory;
import de.gnwi.jdpd.samples.harmonicBonds.ParticlePairHarmonicBond;
import de.gnwi.jdpd.parameters.ParticleTypes;
import de.gnwi.jdpd.utilities.BoxSize;
import de.gnwi.jdpd.utilities.Electrostatics;
import de.gnwi.jdpd.utilities.MoleculeDescription;
import de.gnwi.jdpd.movement.MoleculeFixationDescription;
import de.gnwi.jdpd.rg.MoleculeRgCalculationDescription;
import de.gnwi.jdpd.movement.MoleculeVelocityFixationDescription;
import de.gnwi.jdpd.nearestNeighbor.NearestNeighborBaseParticleDescription;
import de.gnwi.jdpd.utilities.GravitationalAcceleration;
import de.gnwi.jdpd.utilities.PeriodicBoundaries;
import java.util.HashMap;

/**
 * Interface for simulation input
 * 
 * @author Achim Zielesny
 */
public interface IInput {

    // <editor-fold defaultstate="collapsed" desc="Factory (get)">
    /**
     * Returns factory
     * 
     * @return Factory
     */
    Factory getFactory();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Particle description properties (get)">
    /**
     * Particle types
     * 
     * @return Particle types
     */
    ParticleTypes getParticleTypes();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Chemical system description properties (get)">
    /**
     * Molecule descriptions
     * 
     * @return Molecule descriptions
     */
    MoleculeDescription[] getMoleculeDescriptions();
    
    /**
     * Box size in DPD units
     * 
     * @return Box size in DPD units
     */
    BoxSize getBoxSize();

    /**
     * Molecule fixation descriptions
     * 
     * @return Molecule fixation descriptions
     */
    MoleculeFixationDescription[] getMoleculeFixationDescriptions();

    /**
     * Molecule boundary descriptions
     * 
     * @return Molecule boundary descriptions or null if no molecule is bounded
     */
    MoleculeBoundaryDescription[] getMoleculeBoundaryDescriptions();
    
    /**
     * Molecule velocity fixation descriptions
     * 
     * @return Molecule velocity fixation descriptions or null if no molecule velocity is fixed
     */
    MoleculeVelocityFixationDescription[] getMoleculeVelocityFixationDescriptions();

    /**
     * Molecule acceleration descriptions
     * 
     * @return Molecule acceleration descriptions or null if no molecule is accelerated
     */
    MoleculeAccelerationDescription[] getMoleculeAccelerationDescriptions();

    /**
     * Molecule particle radius of gyration (Rg) calculation descriptions
     * 
     * @return Molecule particle radius of gyration (Rg) calculation descriptions
     */
    MoleculeRgCalculationDescription[] getMoleculeRgCalculationDescriptions();
    
    /**
     * Nearest-neighbor particle descriptions
     * 
     * @return Nearest-neighbor particle descriptions
     */
    NearestNeighborBaseParticleDescription[] getNearestNeighborParticleDescriptions();
    
    /**
     * Nearest-neighbor distance in DPD units
     * 
     * @return Nearest-neighbor distance in DPD units
     */
    double getNearestNeighborDistance();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Interaction description properties (get)">
    /**
     * Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * 
     * @return Temperature
     */
    double getTemperature();

    /**
     * DPD sigma parameter in DPD units
     * 
     * @return DPD sigma parameter
     */
    double getDpdSigma();
    
    /**
     * Conservative force parameters a(ij)
     * 
     * @param aParticleTypes Particle types
     * @return Conservative force parameters a(ij)
     */
    double[][] getAij(ParticleTypes aParticleTypes);
    
    /**
     * Hash map that maps particle-pair key to corresponding particle-pair bond.
     * Note: Particle-pair key has to be created with Utils.getParticlePairKey().
     * 
     * @return Hash map that maps particle-pair key to corresponding particle-pair bond
     */
    HashMap<String, ParticlePairHarmonicBond> getParticlePairBondMap();
    
    /**
     * Flag for use of Gaussian random for random force
     * 
     * @return True : Random DPD force is driven by random variable with 
     *                Gaussian distribution with zero mean and unit variance
     *                (slower)
     *         False: Random DPD force is driven by random variable with uniform
     *                distribution with zero mean and unit variance (faster)
     */
    boolean isGaussianRandomDpdForce();
    
    /**
     * Electrostatics parameters
     * 
     * @return Electrostatics parameters or null if none are available
     */
    Electrostatics getElectrostatics();
    
    /**
     * Gravitational acceleration
     * 
     * @return Gravitational acceleration
     */
    GravitationalAcceleration getGravitationalAcceleration();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Simulation description properties (get)">
    /**
     * Time step number for central simulation loop
     * 
     * @return Time step number
     */
    int getTimeStepNumber();

    /**
     * Time step length in DPD units
     * 
     * @return Time step length
     */
    double getTimeStepLength();

    /**
     * Time step frequency for output
     * 
     * @return Time step frequency for output
     */
    int getTimeStepFrequencyForOutput();
    
    /**
     * Number of initial potential minimisation steps
     * 
     * @return Number of initial potential minimisation steps
     */
    int getInitialPotentialEnergyMinimizationStepNumber();
    
    /**
     * Flag for initial potential energy minimisation step output
     * 
     * @return True: Potential energy minimisation step output is generated, false: Otherwise (NO output)
     */
    boolean isInitialPotentialEnergyMinimizationStepOutput();
    
    /**
     * Periodic boundaries
     * 
     * @return Periodic boundaries
     */
    PeriodicBoundaries getPeriodicBoundaries();

    /**
     * Flag for use of DPD unit masses:
     * True : DPD masses of all particles are set to 1
     * False: The DPD mass of the most lightweight particle (often water) is set to 1. 
     *        The masses of all other particles are set in relation to their 
     *        molar mass ratios to the most lightweight particle.
     * 
     * @return Flag value for use of DPD unit masses
     */
    boolean isDpdUnitMass();

    /**
     * True: Velocity scaling is performed for every simulation step, false: Otherwise
     * 
     * @return True: Velocity scaling is performed for every simulation step, false: Otherwise
     */
    boolean isVelocityScaling();
    
    /**
     * Seed for random number generation
     * 
     * @return Seed for random number generation
     */
    int getRandomSeed();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Simulation counts properties (get)">
    /**
     * Particle number
     * 
     * @return Particle number
     */
    int getParticleNumber();
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Additional properties (get)">
    /**
     * Maximum number of position correction trials
     * 
     * @return Maximum number of position correction trials
     */
    int getMaximumNumberOfPositionCorrectionTrials();
    // </editor-fold>
    
}
