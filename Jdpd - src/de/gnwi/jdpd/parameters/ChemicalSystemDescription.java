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
package de.gnwi.jdpd.parameters;

import de.gnwi.jdpd.movement.MoleculeAccelerationInfo;
import de.gnwi.jdpd.movement.MoleculeBoundaryInfo;
import de.gnwi.jdpd.movement.MoleculeVelocityFixationInfo;
import de.gnwi.jdpd.movement.MoleculeFixationInfo;
import de.gnwi.jdpd.movement.MoleculeSphereInfo;
import de.gnwi.jdpd.nearestNeighbor.NearestNeighborManager;
import de.gnwi.jdpd.rg.RgCalculator;
import de.gnwi.jdpd.utilities.BoxSize;

/**
 * Chemical system description parameters
 * 
 * @author Achim Zielesny
 */
public class ChemicalSystemDescription {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Molecule types
     */
    private final MoleculeTypes moleculeTypes;

    /**
     * Box size in DPD units
     */
    private final BoxSize boxSize;
    
    /**
     * Molecule fixation infos
     */
    private final MoleculeFixationInfo[] moleculeFixationInfos;
    
    /**
     * Molecule boundary infos
     */
    private final MoleculeBoundaryInfo[] moleculeBoundaryInfos;
    
    /**
     * Molecule sphere infos
     */
    private final MoleculeSphereInfo[] moleculeSphereInfos;
    
    /**
     * Molecule velocity fixation infos
     */
    private final MoleculeVelocityFixationInfo[] moleculeVelocityFixationInfos;
    
    /**
     * Molecule acceleration infos
     */
    private final MoleculeAccelerationInfo[] moleculeAccelerationInfos;
    
    /**
     * Rg calculators
     */
    private final RgCalculator[] rgCalculators;
    
    /**
     * True: Rg calculation; false: Otherwise
     */
    private final boolean isRgCalculation;
    
    /**
     * Nearest-neighbor manager
     */
    private final NearestNeighborManager nearestNeighborManager;
    
    /**
     * True: Nearest-neighbor particle determination; false: Otherwise
     */
    private final boolean isNearestNeighborParticleDetermination;
    
    /**
     * Nearest-neighbor distance in DPD units
     */
    private final double nearestNeighborDistance;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aMoleculeTypes Molecule types
     * @param aBoxSize Box size in DPD units
     * @param aMoleculeFixationInfos Molecule fixation infos (may be null)
     * @param aMoleculeBoundaryInfos Molecule boundary infos (may be null)
     * @param aMoleculeSphereInfos Molecule sphere infos (may be null)
     * @param aMoleculeVelocityFixationInfos Molecule velocity fixation infos (may be null)
     * @param aMoleculeAccelerationInfos Molecule acceleration infos (may be null)
     * @param aRgCalculators Rg calculators (may be null)
     * @param aNearestNeighborManager (may be null)
     * @param aNearestNeighborDistance Nearest-neighbor distance in DPD units
     */
    public ChemicalSystemDescription(
        MoleculeTypes aMoleculeTypes,
        BoxSize aBoxSize,
        MoleculeFixationInfo[] aMoleculeFixationInfos,
        MoleculeBoundaryInfo[] aMoleculeBoundaryInfos,
        MoleculeSphereInfo[] aMoleculeSphereInfos,
        MoleculeVelocityFixationInfo[] aMoleculeVelocityFixationInfos,
        MoleculeAccelerationInfo[] aMoleculeAccelerationInfos,
        RgCalculator[] aRgCalculators,
        NearestNeighborManager aNearestNeighborManager,
        double aNearestNeighborDistance
    ) {
        this.moleculeTypes = aMoleculeTypes;
        this.boxSize = aBoxSize;
        this.moleculeFixationInfos = aMoleculeFixationInfos;
        this.moleculeBoundaryInfos = aMoleculeBoundaryInfos;
        this.moleculeSphereInfos = aMoleculeSphereInfos;
        this.moleculeVelocityFixationInfos = aMoleculeVelocityFixationInfos;
        this.moleculeAccelerationInfos = aMoleculeAccelerationInfos;
        
        this.rgCalculators = aRgCalculators;
        this.isRgCalculation = this.rgCalculators != null;
        
        this.nearestNeighborManager = aNearestNeighborManager;
        this.isNearestNeighborParticleDetermination = this.nearestNeighborManager != null;
        this.nearestNeighborDistance = aNearestNeighborDistance;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Molecule types
     * 
     * @return Molecule types
     */
    public MoleculeTypes getMoleculeTypes() {
        return this.moleculeTypes;
    }
    
    /**
     * Box size in DPD units
     * 
     * @return Box size in DPD units
     */
    public BoxSize getBoxSize() {
        return this.boxSize;
    }
    
    /**
     * Molecule fixation infos
     * 
     * @return Molecule fixation infos
     */
    public MoleculeFixationInfo[] getMoleculeFixationInfos() {
        return this.moleculeFixationInfos;
    }
    
    /**
     * Molecule boundary infos
     * 
     * @return Molecule boundary infos
     */
    public MoleculeBoundaryInfo[] getMoleculeBoundaryInfos() {
        return this.moleculeBoundaryInfos;
    }
    
    /**
     * Molecule sphere infos
     * 
     * @return Molecule sphere infos
     */
    public MoleculeSphereInfo[] getMoleculeSphereInfos() {
        return this.moleculeSphereInfos;
    }
    
    /**
     * Molecule velocity fixation infos
     * 
     * @return Molecule velocity fixation infos
     */
    public MoleculeVelocityFixationInfo[] getMoleculeVelocityFixationInfos() {
        return this.moleculeVelocityFixationInfos;
    }
    
    /**
     * Molecule acceleration infos
     * 
     * @return Molecule acceleration infos
     */
    public MoleculeAccelerationInfo[] getMoleculeAccelerationInfos() {
        return this.moleculeAccelerationInfos;
    }
    
    /**
     * Rg calculators
     * 
     * @return Rg calculators
     */
    public RgCalculator[] getRgCalculators() {
        return this.rgCalculators;
    }
    
    /**
     * True: Rg calculation; false: Otherwise
     * 
     * @return True: Rg calculation; false: Otherwise
     */
    public boolean isRgCalculation() {
        return this.isRgCalculation;
    }
        
    /**
     * Nearest-neighbor manager
     * 
     * @return Nearest-neighbor manager or null if none exists
     */
    public NearestNeighborManager getNearestNeighborManager() {
        return this.nearestNeighborManager;
    }

    /**
     * True: Nearest-neighbor particle determination; false: Otherwise
     * 
     * @return True: Nearest-neighbor particle determination; false: Otherwise
     */
    public boolean isNearestNeighborParticleDetermination() {
        return this.isNearestNeighborParticleDetermination;
    }
    
    /**
     * Nearest-neighbor distance in DPD units
     * 
     * @return Nearest-neighbor distance in DPD units
     */
    public double getNearestNeighborDistance() {
        return this.nearestNeighborDistance;
    }
    // </editor-fold>
    
}
