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

import de.gnwi.jdpdsp.samples.harmonicBonds.HarmonicBondChunkArrays;
import java.util.LinkedList;

/**
 * Particle arrays
 * 
 * @author Achim Zielesny
 */
public class ParticleArrays {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Current x-components of particle positions in simulation box
     */
    private final float[] r_x;

    /**
     * Current y-components of particle positions in simulation box
     */
    private final float[] r_y;
    
    /**
     * Current z-components of particle positions in simulation box
     */
    private final float[] r_z;

    /**
     * Old x-components of particle positions in simulation box
     */
    private final float[] rOld_x;

    /**
     * Old y-components of particle positions in simulation box
     */
    private final float[] rOld_y;
    
    /**
     * Old z-components of particle positions in simulation box
     */
    private final float[] rOld_z;

    /**
     * Current x-components of particle velocities
     */
    private final float[] v_x;

    /**
     * Current y-components of particle velocities
     */
    private final float[] v_y;
    
    /**
     * Current z-components of particle velocities
     */
    private final float[] v_z;

    /**
     * New x-components of particle velocities
     */
    private final float[] vNew_x;

    /**
     * New y-components of particle velocities
     */
    private final float[] vNew_y;
    
    /**
     * New z-components of particle velocities
     */
    private final float[] vNew_z;

    /**
     * Current x-components of particle forces
     */
    private final float[] f_x;

    /**
     * Current y-components of particle forces
     */
    private final float[] f_y;
    
    /**
     * Current z-components of particle forces
     */
    private final float[] f_z;

    /**
     * Force two x-component
     */
    private final float[] fTwo_x;

    /**
     * Force two y-component
     */
    private final float[] fTwo_y;
    
    /**
     * Force two z-component
     */
    private final float[] fTwo_z;
    
    /**
     * Tokens of particles (e.g. "H2O")
     */
    private final String[] particleTokens;
    
    /**
     * Molecule names of particles
     */
    private final String[] moleculeNamesOfParticles;
    
    /**
     * Indices of particle type
     */
    private final int[] particleTypeIndices;
    
    /**
     * Indices of molecule type
     */
    private final int[] moleculeTypeIndices;
    
    /**
     * Indices of molecules
     */
    private final int[] moleculeIndices;
    
    /**
     * Charges of particles
     */
    private final float[] charges;
    
    /**
     * DPD masses of particles
     */
    private final float[] dpdMasses;
    
    /**
     * Molar masses of particles
     */
    private final float[] molarMasses;

    /**
     * List of bond chunk arrays
     */
    private final LinkedList<HarmonicBondChunkArrays> bondChunkArraysList;
    
    /**
     * Charged particle indices
     */
    private final int[] chargedParticleIndices;

    /**
     * Charged particles x-positions
     */
    private final float[] chargedParticles_r_x;

    /**
     * Charged particles y-positions
     */
    private final float[] chargedParticles_r_y;

    /**
     * Charged particles z-positions
     */
    private final float[] chargedParticles_r_z;

    /**
     * True: Nearest-neighbor for base particle is to be determined, false: Otherwise
     */
    private final boolean[] isNearestNeighborBaseParticleDeterminations;
    
    /**
     * Distances to nearest (non-bonded) neighbor particles (of another molecule)
     */
    private final float[] nearestNeighborDistances;

    /**
     * Indices of nearest (non-bonded) neighbor particles (of another molecule)
     */
    private final int[] nearestNeighborParticleIndices;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aR_x Current x-components of particle positions in simulation box
     * @param aR_y Current y-components of particle positions in simulation box
     * @param aR_z Current z-components of particle positions in simulation box
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aParticleTokens Tokens of particles (e.g. "H2O")
     * @param aMoleculeNamesOfParticles Molecule names of particles
     * @param aParticleTypeIndices Indices of particle type
     * @param aMoleculeTypeIndices Indices of molecule type
     * @param aMoleculeIndices Indices of molecules
     * @param aCharges Charges of particles
     * @param aDpdMasses DPD masses of particles
     * @param aMolarMasses Molar masses of particles
     * @param aBondChunkArraysList List of bond chunk arrays (may be null)
     * @param aChargedParticleIndices Charged particle indices (may be null)
     * @param anIsNearestNeighborBaseParticleDeterminations True: Nearest-neighbor 
     * for base particle is to be determined, false: Otherwise (may be null)
     * @param aNearestNeighborDistances Distances to nearest (non-bonded) 
     * neighbor particles (of another molecule) (may be null)
     * @param aNearestNeighborParticleIndices Indices of nearest (non-bonded) 
     * neighbor particles (of another molecule) (may be null)
     */
    public ParticleArrays(
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        String[] aParticleTokens,
        String[] aMoleculeNamesOfParticles,
        int[] aParticleTypeIndices,
        int[] aMoleculeTypeIndices,
        int[] aMoleculeIndices,
        float[] aCharges,
        float[] aDpdMasses,
        float[] aMolarMasses,
        LinkedList<HarmonicBondChunkArrays> aBondChunkArraysList,
        int[] aChargedParticleIndices,
        boolean[] anIsNearestNeighborBaseParticleDeterminations,
        float[] aNearestNeighborDistances,
        int[] aNearestNeighborParticleIndices
    ) {
        this.r_x = aR_x;
        this.r_y = aR_y;
        this.r_z = aR_z;

        this.v_x = aV_x;
        this.v_y = aV_y;
        this.v_z = aV_z;

        this.rOld_x = new float[this.r_x.length];
        this.rOld_y = new float[this.r_x.length];
        this.rOld_z = new float[this.r_x.length];

        this.vNew_x = new float[this.r_x.length];
        this.vNew_y = new float[this.r_x.length];
        this.vNew_z = new float[this.r_x.length];
        
        this.f_x = new float[this.r_x.length];
        this.f_y = new float[this.r_x.length];
        this.f_z = new float[this.r_x.length];

        this.fTwo_x = new float[this.r_x.length];
        this.fTwo_y = new float[this.r_x.length];
        this.fTwo_z = new float[this.r_x.length];
        
        this.particleTokens = aParticleTokens;
        this.moleculeNamesOfParticles = aMoleculeNamesOfParticles;
        this.particleTypeIndices = aParticleTypeIndices;
        this.moleculeTypeIndices = aMoleculeTypeIndices;
        this.moleculeIndices = aMoleculeIndices;
        this.charges = aCharges;
        this.dpdMasses = aDpdMasses;
        this.molarMasses = aMolarMasses;
        this.bondChunkArraysList = aBondChunkArraysList;
        this.chargedParticleIndices = aChargedParticleIndices;
        
        if (this.chargedParticleIndices != null) {
            this.chargedParticles_r_x = new float[this.chargedParticleIndices.length];
            this.chargedParticles_r_y = new float[this.chargedParticleIndices.length];
            this.chargedParticles_r_z = new float[this.chargedParticleIndices.length];
        } else {
            this.chargedParticles_r_x = null;
            this.chargedParticles_r_y = null;
            this.chargedParticles_r_z = null;
        }
        
        this.isNearestNeighborBaseParticleDeterminations = anIsNearestNeighborBaseParticleDeterminations;
        this.nearestNeighborDistances = aNearestNeighborDistances;
        this.nearestNeighborParticleIndices = aNearestNeighborParticleIndices;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Returns if charged particles exist
     * 
     * @return True: Charged particles exist, false: Otherwise
     */
    public boolean hasChargedParticles() {
        return this.chargedParticleIndices != null;
    }
    
    /**
     * Returns charged particles x-positions
     * 
     * @return Charged particles x-positions or null if none exist
     */
    public float[] getChargedParticles_r_x() {
        if (this.hasChargedParticles()) {
            for (int i = 0; i < this.chargedParticleIndices.length; i++) {
                this.chargedParticles_r_x[i] = this.r_x[this.chargedParticleIndices[i]];
            }
            return this.chargedParticles_r_x;
        } else {
            return null;
        }
    }
    
    /**
     * Returns charged particles y-positions
     * 
     * @return Charged particles y-positions or null if none exist
     */
    public float[] getChargedParticles_r_y() {
        if (this.hasChargedParticles()) {
            for (int i = 0; i < this.chargedParticleIndices.length; i++) {
                this.chargedParticles_r_y[i] = this.r_y[this.chargedParticleIndices[i]];
            }
            return this.chargedParticles_r_y;
        } else {
            return null;
        }
    }
    
    /**
     * Returns charged particles z-positions
     * 
     * @return Charged particles z-positions or null if none exist
     */
    public float[] getChargedParticles_r_z() {
        if (this.hasChargedParticles()) {
            for (int i = 0; i < this.chargedParticleIndices.length; i++) {
                this.chargedParticles_r_z[i] = this.r_z[this.chargedParticleIndices[i]];
            }
            return this.chargedParticles_r_z;
        } else {
            return null;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Current x-components of particle positions in simulation box
     * 
     * @return Current x-components of particle positions in simulation box
     */
    public float[] getR_x() {
        return this.r_x;
    }

    /**
     * Current y-components of particle positions in simulation box
     * 
     * @return Current y-components of particle positions in simulation box
     */
    public float[] getR_y() {
        return this.r_y;
    }

    /**
     * Current z-components of particle positions in simulation box
     * 
     * @return Current z-components of particle positions in simulation box
     */
    public float[] getR_z() {
        return this.r_z;
    }

    /**
     * Old x-components of particle positions in simulation box
     * 
     * @return Old x-components of particle positions in simulation box
     */
    public float[] getRold_x() {
        return this.rOld_x;
    }

    /**
     * Old y-components of particle positions in simulation box
     * 
     * @return Old y-components of particle positions in simulation box
     */
    public float[] getRold_y() {
        return this.rOld_y;
    }

    /**
     * Old z-components of particle positions in simulation box
     * 
     * @return Old z-components of particle positions in simulation box
     */
    public float[] getRold_z() {
        return this.rOld_z;
    }
    
    /**
     * Current x-components of particle velocities
     * 
     * @return Current x-components of particle velocities
     */
    public float[] getV_x() {
        return this.v_x;
    }

    /**
     * Current y-components of particle velocities
     * 
     * @return Current y-components of particle velocities
     */
    public float[] getV_y() {
        return this.v_y;
    }

    /**
     * Current z-components of particle velocities
     * 
     * @return Current z-components of particle velocities
     */
    public float[] getV_z() {
        return this.v_z;
    }
    
    /**
     * New x-components of particle velocities
     * 
     * @return New x-components of particle velocities
     */
    public float[] getVnew_x() {
        return this.vNew_x;
    }
    
    /**
     * New y-components of particle velocities
     * 
     * @return New y-components of particle velocities
     */
    public float[] getVnew_y() {
        return this.vNew_y;
    }
    
    /**
     * New z-components of particle velocities
     * 
     * @return New z-components of particle velocities
     */
    public float[] getVnew_z() {
        return this.vNew_z;
    }
    
    /**
     * Current x-components of particle forces
     * 
     * @return Current x-components of particle forces
     */
    public float[] getF_x() {
        return this.f_x;
    }
    
    /**
     * Current y-components of particle forces
     * 
     * @return Current y-components of particle forces
     */
    public float[] getF_y() {
        return this.f_y;
    }
    
    /**
     * Current z-components of particle forces
     * 
     * @return Current z-components of particle forces
     */
    public float[] getF_z() {
        return this.f_z;
    }
    
    /**
     * Force two x-component
     * 
     * @return Force two x-component
     */
    public float[] getFtwo_x() {
        return this.fTwo_x;
    }
    
    /**
     * Force two y-component
     * 
     * @return Force two y-component
     */
    public float[] getFtwo_y() {
        return this.fTwo_y;
    }
    
    /**
     * Force two z-component
     * 
     * @return Force two z-component
     */
    public float[] getFtwo_z() {
        return this.fTwo_z;
    }
    
    /**
     * Tokens of particles (e.g. "H2O")
     * 
     * @return Tokens of particles (e.g. "H2O")
     */
    public String[] getParticleTokens() {
        return this.particleTokens;
    }
    
    /**
     * Molecule names of particles
     * 
     * @return Molecule names of particles
     */
    public String[] getMoleculeNamesOfParticles() {
        return this.moleculeNamesOfParticles;
    }
    
    /**
     * Indices of particle type
     * 
     * @return Indices of particle type
     */
    public int[] getParticleTypeIndices() {
        return this.particleTypeIndices;
    }
    
    /**
     * Indices of molecule type
     * 
     * @return Indices of molecule type
     */
    public int[] getMoleculeTypeIndices() {
        return this.moleculeTypeIndices;
    }
    
    /**
     * Indices of molecules
     * 
     * @return Indices of molecules
     */
    public int[] getMoleculeIndices() {
        return this.moleculeIndices;
    }
    
    /**
     * Charges of particles
     * 
     * @return Charges of particles
     */
    public float[] getCharges() {
        return this.charges;
    }
    
    /**
     * DPD masses of particles
     * 
     * @return DPD masses of particles
     */
    public float[] getDpdMasses() {
        return this.dpdMasses;
    }
    
    /**
     * Molar masses of particles
     * 
     * @return Molar masses of particles
     */
    public float[] getMolarMasses() {
        return this.molarMasses;
    }

    /**
     * List of bond chunk arrays
     * 
     * @return List of bond chunk arrays
     */
    public LinkedList<HarmonicBondChunkArrays> getBondChunkArraysList() {
        return this.bondChunkArraysList;
    }
    
    /**
     * Returns if bonds exist
     * 
     * @return True: Bonds exist (and may be accessed with method getBondChunkArraysList()), false: Otherwise
     */
    public boolean hasBonds() {
        return this.bondChunkArraysList != null;
    }

    /**
     * Charged particle indices
     * 
     * @return Charged particle indices
     */
    public int[] getChargedParticleIndices() {
        return this.chargedParticleIndices;
    }

    /**
     * True: Nearest-neighbor for base particle is to be determined, false: Otherwise
     * 
     * @return  True: Nearest-neighbor for base particle is to be determined, 
     * false: Otherwise, or null if none exist
     */
    public boolean[] isNearestNeighborBaseParticleDeterminations() {
        return this.isNearestNeighborBaseParticleDeterminations;
    };
    
    /**
     * Distances to nearest (non-bonded) neighbor particles (of another molecule)
     * 
     * @return Distances to nearest (non-bonded) neighbor particles (of another 
     * molecule) or null if none exist
     */
    public float[] getNearestNeighborDistances() {
        return this.nearestNeighborDistances;
    }

    /**
     * Indices of nearest (non-bonded) neighbor particles (of another molecule)
     * 
     * @return Indices of nearest (non-bonded) neighbor particles (of another 
     * molecule) or null if none exist
     */
    public int[] getNearestNeighborParticleIndices() {
        return this.nearestNeighborParticleIndices;
    }
    // </editor-fold>
    
}
