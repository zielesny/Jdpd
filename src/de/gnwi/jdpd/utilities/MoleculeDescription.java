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
package de.gnwi.jdpd.utilities;

import de.gnwi.jdpd.samples.harmonicBonds.HarmonicBond;

/**
 * Molecule description
 * 
 * @author Achim Zielesny
 */
public class MoleculeDescription {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Name of molecule
     */
    private final String moleculeName;

    /**
     * Total number of all particles of molecules
     */
    private final int totalMoleculeParticleNumber;

    /**
     * Number of particles of a single molecule
     */
    private final int singleMoleculeParticleNumber;

    /**
     * Current x-position of particle in simulation box
     */
    private final double[] r_x;

    /**
     * Current y-position of particle in simulation box
     */
    private final double[] r_y;
    
    /**
     * Current z-position of particle in simulation box
     */
    private final double[] r_z;

    /**
     * Particle tokens
     */
    private final String[] particleTokens;
    
    /**
     * Particle backbone indices
     */
    private final int[] particleBackboneIndices;
    
    /**
     * Bond offsets
     */
    private final int[][] bondOffsets;
    
    /**
     * Backbone bonds
     */
    private final HarmonicBond[] backboneBonds;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aMoleculeName Name of molecule
     * @param aTotalMoleculeParticleNumber Total number of all particles of molecules
     * @param aSingleMoleculeParticleNumber Number of particles of a single molecule
     * @param aParticleTokens Particle tokens
     * @param aParticleBackboneIndices Particle backbone indices
     * @param aR_x Current x-position of particle in simulation box
     * @param aR_y Current y-position of particle in simulation box
     * @param aR_z Current z-position of particle in simulation box
     * @param aBondOffsets Bond offsets
     * @param aBackboneBonds Backbone bonds (may be null)
     */
    public MoleculeDescription(
        String aMoleculeName,
        int aTotalMoleculeParticleNumber,
        int aSingleMoleculeParticleNumber,
        String[] aParticleTokens,
        int[] aParticleBackboneIndices,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        int[][] aBondOffsets,
        HarmonicBond[] aBackboneBonds) {
        this.moleculeName = aMoleculeName;
        this.totalMoleculeParticleNumber = aTotalMoleculeParticleNumber;
        this.singleMoleculeParticleNumber = aSingleMoleculeParticleNumber;
        this.particleTokens = aParticleTokens;
        this.particleBackboneIndices = aParticleBackboneIndices;
        this.r_x = aR_x;
        this.r_y = aR_y;
        this.r_z = aR_z;
        this.bondOffsets = aBondOffsets;
        this.backboneBonds = aBackboneBonds;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Name of molecule
     * 
     * @return Name of molecule
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }

    /**
     * Total number of all particles of molecules
     * 
     * @return Total number of all particles of molecules
     */
    public int getTotalMoleculeParticleNumber() {
        return this.totalMoleculeParticleNumber;
    }

    /**
     * Number of particles of a single molecule
     * 
     * @return Number of particles of a single molecule
     */
    public int getSingleMoleculeParticleNumber() {
        return this.singleMoleculeParticleNumber;
    }

    /**
     * Particle tokens
     * 
     * @return Particle tokens
     */
    public String[] getParticleTokens() {
        return this.particleTokens;
    };

    /**
     * Particle backbone indices
     * 
     * @return Particle backbone indices
     */
    public int[] getParticleBackboneIndices() {
        return this.particleBackboneIndices;
    };

    /**
     * Current x-position of particle in simulation box
     * 
     * @return Current x-position of particle in simulation box
     */
    public double[] getR_x() {
        return this.r_x;
    };

    /**
     * Current y-position of particle in simulation box
     * 
     * @return Current y-position of particle in simulation box
     */
    public double[] getR_y() {
        return this.r_y;
    };
    
    /**
     * Current z-position of particle in simulation box
     * 
     * @return Current z-position of particle in simulation box
     */
    public double[] getR_z() {
        return this.r_z;
    };
    
    /**
     * Bond offsets
     * 
     * @return Bond offsets
     */
    public int[][] getBondOffsets() {
        return this.bondOffsets;
    }

    /**
     * Backbone bonds
     * 
     * @return Backbone bonds
     */
    public HarmonicBond[] getBackboneBonds() {
        return this.backboneBonds;
    }

    /**
     * Number of backbone bonds
     * 
     * @return Number of backbone bonds
     */
    public int getBackboneBondNumber() {
        if (this.backboneBonds == null) {
            return 0;
        } else {
            return this.backboneBonds.length;
        }
    }
    
    /**
     * Number of molecules
     * 
     * @return Number of molecules
     */
    public int getMoleculeNumber() {
        return this.totalMoleculeParticleNumber/this.singleMoleculeParticleNumber;
    }
    // </editor-fold>
    
}
