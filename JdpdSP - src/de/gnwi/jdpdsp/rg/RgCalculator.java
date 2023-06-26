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
package de.gnwi.jdpdsp.rg;

import de.gnwi.jdpdsp.utilities.BoxSize;
import de.gnwi.jdpdsp.utilities.PeriodicBoundaries;
import de.gnwi.jdpdsp.utilities.Utils;
import org.apache.commons.math3.util.FastMath;

/**
 * Particle radius of gyration calculation
 * 
 * @author Achim Zielesny
 */
public class RgCalculator {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Molecule name
     */
    private final String moleculeName;

    /**
     * Start index of molecule type in aR_x etc. (see method getMeanRadiusOfGyration())
     */
    private final int moleculeTypeStartIndex;
    
    /**
     * Molar mass
     */
    private final float molarMass;
    
    /**
     * Number of particles per molecule
     */
    private final int singleMoleculeParticleNumber;
    
    /**
     * Total number of particles for all molecules
     */
    private final int totalMoleculeParticleNumber;
    
    /**
     * Box size
     */
    private final BoxSize boxSize;
    
    /**
     * Periodic boundaries
     */
    private final PeriodicBoundaries periodicBoundaries;
    
    /**
     * Number of molecules
     */
    private final int numberOfMolecules;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * X-components of particle positions in simulation box
     */
    private float[] r_x;

    /**
     * Y-components of particle positions in simulation box
     */
    private float[] r_y;
    
    /**
     * Z-components of particle positions in simulation box
     */
    private float[] r_z;
    
    /**
     * Molar masses of particles
     * NOTE: Array corresponds to r_x, r_y, r_z
     */
    private float[] particleMolarMasses;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aMoleculeTypeStartIndex Start index of molecule type in aR_x etc. (see method getMeanRadiusOfGyration())
     * @param aMoleculeName Molecule name
     * @param aMolarMass Molar mass of molecule
     * @param aSingleMoleculeParticleNumber Number of particles per molecule
     * @param aTotalMoleculeParticleNumber Total number of particles for all molecules
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public RgCalculator(
        int aMoleculeTypeStartIndex,
        String aMoleculeName,
        float aMolarMass,
        int aSingleMoleculeParticleNumber,
        int aTotalMoleculeParticleNumber,
        BoxSize aBoxSize,
        PeriodicBoundaries aPeriodicBoundaries
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aMoleculeTypeStartIndex < 0) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aMoleculeTypeStartIndex is illegal.");
        }
        if (aMoleculeName == null || aMoleculeName.isEmpty()) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aMoleculeName is null/empty.");
        }
        if (aMolarMass <= 0.0f) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aMolarMass is illegal.");
        }
        if (aSingleMoleculeParticleNumber < 1) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aSingleMoleculeParticleNumber is illegal.");
        }
        if (aTotalMoleculeParticleNumber < 1) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aTotalMoleculeParticleNumber is illegal.");
        }
        if (aTotalMoleculeParticleNumber % aSingleMoleculeParticleNumber != 0) {
            throw new IllegalArgumentException("RgCalculator.Constructor: Ratio aTotalMoleculeParticleNumber/aSingleMoleculeParticleNumber is illegal.");
        }
        if (aBoxSize == null) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aBoxSize is null.");
        }
        if (aPeriodicBoundaries == null) {
            throw new IllegalArgumentException("RgCalculator.Constructor: aPeriodicBoundaries is null.");
        }
        // </editor-fold>
        this.moleculeTypeStartIndex = aMoleculeTypeStartIndex;
        this.moleculeName = aMoleculeName;
        this.molarMass = aMolarMass;
        this.singleMoleculeParticleNumber = aSingleMoleculeParticleNumber;
        this.totalMoleculeParticleNumber = aTotalMoleculeParticleNumber;
        this.boxSize = aBoxSize;
        this.periodicBoundaries = aPeriodicBoundaries;
        this.numberOfMolecules = this.totalMoleculeParticleNumber/this.singleMoleculeParticleNumber;
        this.r_x = null;
        this.r_y = null;
        this.r_z = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Calculates the particle radius of gyration
     * (No checks are performed)
     * NOTE: aR_x, aR_y, aR_z and aParticleMolarMasses correspond to each other,
     *       e.g. aR_x[i] and aParticleMolarMasses[i] refer to the same 
     *       particle i
     * 
     * @param aR_x x positions of particles in simulation box
     * @param aR_y y positions of particles in simulation box
     * @param aR_z z positions of particles in simulation box
     * @param aParticleMolarMasses Molar masses of particles in simulation box
     * @return Mean particle radius of gyration of molecules
     */
    public float getMeanRadiusOfGyration(
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aParticleMolarMasses
    ) {
        if (this.r_x == null) {
            this.r_x = new float[aR_x.length];
            this.r_y = new float[aR_y.length];
            this.r_z = new float[aR_z.length];
        }
        // IMPORTANT: Copy aR_x, aR_y, aR_z to this.r_x, this.r_y, this.r_z
        // since the latter arrays may be changed during Rg calculation
        Utils.copyToOld(
            aR_x, 
            aR_y, 
            aR_z, 
            this.r_x, 
            this.r_y, 
            this.r_z
        );
        this.particleMolarMasses = aParticleMolarMasses;
        float tmpRgSum = 0;
        for (int i = 0; i < this.numberOfMolecules; i++) {
            int tmpSingleMoleculeStartIndex = this.moleculeTypeStartIndex + i * this.singleMoleculeParticleNumber;
            int tmpSingleMoleculeExclusiveEndIndex = this.moleculeTypeStartIndex + (i + 1) * this.singleMoleculeParticleNumber;
            tmpRgSum += this.getRadiusOfGyrationOfSingleMolecule(tmpSingleMoleculeStartIndex, tmpSingleMoleculeExclusiveEndIndex);
        }
        return tmpRgSum / (float) this.numberOfMolecules;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Molecule name
     * 
     * @return Molecule name
     */
    public String getMoleculeName() {
        return this.moleculeName;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Calculates the particle radius of gyration of a single molecule with
     * specified start and (exclusive) end index in this.r_x etc.
     * (No checks are performed)
     * 
     * @param aSingleMoleculeStartIndex Single molecule start index in this.r_x etc.
     * @param aSingleMoleculeExclusiveEndIndex Exclusive single molecule end index in this.r_x etc.
     * @return Particle radius of gyration of single molecule
     */
    private float getRadiusOfGyrationOfSingleMolecule(int aSingleMoleculeStartIndex, int aSingleMoleculeExclusiveEndIndex) {
        // <editor-fold defaultstate="collapsed" desc="Correct for periodic boundaries">
        if (this.periodicBoundaries.isPeriodicBoundary()) {
            float tmpXmin = Float.MAX_VALUE;
            float tmpXmax = Float.MIN_VALUE;
            float tmpYmin = Float.MAX_VALUE;
            float tmpYmax = Float.MIN_VALUE;
            float tmpZmin = Float.MAX_VALUE;
            float tmpZmax = Float.MIN_VALUE;
            for (int i = aSingleMoleculeStartIndex; i < aSingleMoleculeExclusiveEndIndex; i++) {
                if (this.r_x[i] < tmpXmin) {
                    tmpXmin = this.r_x[i];
                }
                if (this.r_x[i] > tmpXmax) {
                    tmpXmax = this.r_x[i];
                }
                if (this.r_y[i] < tmpYmin) {
                    tmpYmin = this.r_y[i];
                }
                if (this.r_y[i] > tmpYmax) {
                    tmpYmax = this.r_y[i];
                }
                if (this.r_z[i] < tmpZmin) {
                    tmpZmin = this.r_z[i];
                }
                if (this.r_z[i] > tmpZmax) {
                    tmpZmax = this.r_z[i];
                }
            }
            float tmpDeltaX = tmpXmax - tmpXmin;
            if (tmpDeltaX > this.boxSize.getXHalfLength()) {
                this.correctR_i(aSingleMoleculeStartIndex, aSingleMoleculeExclusiveEndIndex, this.r_x, this.boxSize.getXHalfLength(), this.boxSize.getXLength());
            }
            float tmpDeltaY = tmpYmax - tmpYmin;
            if (tmpDeltaY > this.boxSize.getYHalfLength()) {
                this.correctR_i(aSingleMoleculeStartIndex, aSingleMoleculeExclusiveEndIndex, this.r_y, this.boxSize.getYHalfLength(), this.boxSize.getYLength());
            }
            float tmpDeltaZ = tmpZmax - tmpZmin;
            if (tmpDeltaZ > this.boxSize.getZHalfLength()) {
                this.correctR_i(aSingleMoleculeStartIndex, aSingleMoleculeExclusiveEndIndex, this.r_z, this.boxSize.getZHalfLength(), this.boxSize.getZLength());
            }
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Calculate center of mass">
        float tmpR_CenterOfMassX = 0.0f;
        float tmpR_CenterOfMassY = 0.0f;
        float tmpR_CenterOfMassZ = 0.0f;
        for (int i = aSingleMoleculeStartIndex; i < aSingleMoleculeExclusiveEndIndex; i++) {
            tmpR_CenterOfMassX += this.r_x[i] * this.particleMolarMasses[i];
            tmpR_CenterOfMassY += this.r_y[i] * this.particleMolarMasses[i];
            tmpR_CenterOfMassZ += this.r_z[i] * this.particleMolarMasses[i];
        }
        tmpR_CenterOfMassX /= this.molarMass;
        tmpR_CenterOfMassY /= this.molarMass;
        tmpR_CenterOfMassZ /= this.molarMass;
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Calculate radius of gyration">
        float tmpRg = 0.0f;
        for (int i = aSingleMoleculeStartIndex; i < aSingleMoleculeExclusiveEndIndex; i++) {
            float tmpDeltaX = this.r_x[i] - tmpR_CenterOfMassX;
            float tmpDeltaY = this.r_y[i] - tmpR_CenterOfMassY;
            float tmpDeltaZ = this.r_z[i] - tmpR_CenterOfMassZ;
            float tmpDelta = tmpDeltaX * tmpDeltaX + tmpDeltaY * tmpDeltaY + tmpDeltaZ * tmpDeltaZ;
            tmpRg += tmpDelta * this.particleMolarMasses[i];
        }
        // </editor-fold>
        return (float) FastMath.sqrt(tmpRg / this.molarMass);
    }
    
    /**
     * Corrects this.r_x, this.r_y or this.r_z
     * (No checks are performed)
     * 
     * @param aStartIndex Start index in aR_i
     * @param anExclusiveEndIndex Exclusive end index in aR_i
     * @param aR_i this.r_x etc.
     * @param aBoxSizeHalfLength this.boxSize.getXHalfLength() etc.
     * @param aBoxSizeLength this.boxSize.getXLength() etc.
     */
    private void correctR_i(
        int aStartIndex, 
        int anExclusiveEndIndex,
        float[] aR_i,
        float aBoxSizeHalfLength,
        float aBoxSizeLength
    ) {
        float tmpRiGreaterHalfLength = Float.MAX_VALUE;
        for (int i = aStartIndex; i < anExclusiveEndIndex; i++) {
            if (aR_i[i] > aBoxSizeHalfLength) {
                if (aR_i[i] < tmpRiGreaterHalfLength) {
                    tmpRiGreaterHalfLength = aR_i[i];
                }
            }
        }
        for (int i = aStartIndex; i < anExclusiveEndIndex; i++) {
            aR_i[i] -= tmpRiGreaterHalfLength;
            if (aR_i[i] < 0.0f) {
                aR_i[i] += aBoxSizeLength;
            }
        }
    }
    // </editor-fold>

}
