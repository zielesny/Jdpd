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
package de.gnwi.jdpd.utilities;

import de.gnwi.jdpd.interfaces.ILogger;
import java.util.Arrays;

/**
 * Parallelisation-safe cell-partitioned box as a basis for pairwise interaction 
 * calculations with finite cut-off length.
 * 
 * For basic concept and discussion see
 * M.P. Allen and D.J. Tildesley
 * Computer Simulation of Liquids
 * Clarendon Press
 * Oxford 1987
 * Chapter 5.3
 * 
 * @author Achim Zielesny
 */
public class CellBox {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    private static final double TWO = 2.0;
    private static final int INT_TWO = 2;
    private static final double THREE = 3.0;
    private static final int INT_THREE = 3;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * cellNumberX * cellNumberY;
     */
    private final int cellNumberXtimesY;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected final class variables">
    /**
     * Simulation simulationLogger
     */
    protected final ILogger simulationLogger;
    
    /**
     * Box size
     */
    protected final BoxSize boxSize;
    
    /**
     * Periodic boundaries
     */
    protected final PeriodicBoundaries periodicBoundaries;

    /**
     * Cut-off length for partitioning of the box
     */
    protected final double cutOffLength;

    /**
     * Cut-off length for partitioning of the box in x-direction
     */
    protected final double cutOffLengthX;

    /**
     * Cut-off length for partitioning of the box in y-direction
     */
    protected final double cutOffLengthY;

    /**
     * Cut-off length for partitioning of the box in z-direction
     */
    protected final double cutOffLengthZ;
    
    /**
     * Number of cells in x-direction
     */
    protected final int cellNumberX;
    
    /**
     * Number of cells in y-direction
     */
    protected final int cellNumberY;
    
    /**
     * Number of cells in z-direction
     */
    protected final int cellNumberZ;
    
    /**
     * Number of cells that partition box
     */
    protected final int cellNumber;
    
    /**
     * True: Cell-based calculation may be performed, 
     * false: Full loop over particle index pairs is performed
     */
    protected final boolean isCellBasedCalculation;

    /**
     * Cell neighbours
     */
    protected final int[][] cellNeigbours;
    
    /**
     * Cell chunks that are parallelisation-safe
     */
    protected final int[][] parallelizationSafeCellChunks;
    
    /**
     * Size of parallelisation-safe cell chunks
     */
    protected final int cellChunkSize;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructors">
    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation simulationLogger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected CellBox(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength
    ) throws IllegalArgumentException 
    {
        this(
            aSimulationLogger, 
            aBoxSize, 
            aPeriodicBoundaries, 
            aCutOffLength, 
            null, 
            null
        );
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        this.simulationLogger.appendMethodCall("CellBox.Constructor WITHOUT aCellNeigbours, aParallelizationSafeCellChunks");
        // </editor-fold>
    }

    /**
     * Constructor
     * 
     * @param aSimulationLogger Simulation simulationLogger
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries
     * @param aCutOffLength Cut-off length for partitioning of the box
     * @param aCellNeigbours Cell neighbours (may be null then cell neighbours are determined)
     * @param aParallelizationSafeCellChunks Cell chunks that are parallelisation-safe (may be null then cell chunks are determined)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    protected CellBox(
        ILogger aSimulationLogger, 
        BoxSize aBoxSize, 
        PeriodicBoundaries aPeriodicBoundaries, 
        double aCutOffLength,
        int[][] aCellNeigbours,
        int[][] aParallelizationSafeCellChunks
    ) throws IllegalArgumentException 
    {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("CellBox.Constructor: aSimulationLogger is null.");
        }
        if (aBoxSize == null) {
            throw new IllegalArgumentException("CellBox.Constructor: aBoxSize is null.");
        }
        if (aPeriodicBoundaries == null) {
            throw new IllegalArgumentException("CellBox.Constructor: aPeriodicBoundaries is null.");
        }
        if (aCutOffLength <= 0.0) {
            throw new IllegalArgumentException("CellBox.Constructor: aCutOffLength is smaller/equal zero.");
        }
        if (aBoxSize.getXLength() <= aCutOffLength) {
            throw new IllegalArgumentException("CellBox.Constructor: Combination aBoxSize.getXLength() and aForceCutOffLength is illegal.");
        }
        if (aBoxSize.getYLength() <= aCutOffLength) {
            throw new IllegalArgumentException("CellBox.Constructor: Combination aBoxSize.getYLength() and aForceCutOffLength is illegal.");
        }
        if (aBoxSize.getZLength() <= aCutOffLength) {
            throw new IllegalArgumentException("CellBox.Constructor: Combination aBoxSize.getZLength() and aForceCutOffLength is illegal.");
        }
        // </editor-fold>
        this.simulationLogger = aSimulationLogger;
        try {
            this.boxSize = aBoxSize;
            this.periodicBoundaries = aPeriodicBoundaries;
            this.cutOffLength = aCutOffLength;
            // Cut-off lengths in x,y,z-direction that allow for later
            // parallelization:
            // Box lengths divided by these cut-off lengths are a multiple of 3
            // in x and z direction, and a multiple of 2 in y-direction.
            // NOTE: Each cell has 26 neighbors. To avoid double counting of 
            // later pair interactions only 13 cells are specified in method
            // getCellNeighbours(). Thus in y direction a multiple of 2 is 
            // possible in order to get bigger parallelisation safe cell chunks.
            // Multiple of 3 in x direction:
            int tmpCellNumberThirdX = (int) (this.boxSize.getXLength()/(THREE * this.cutOffLength));
            this.cellNumberX = INT_THREE * tmpCellNumberThirdX;
            if (this.cellNumberX%INT_THREE != 0) {
                throw new ArithmeticException("CellBox.initialise: this.cellNumberX can not be divided by 3 without rest");
            }
            this.cutOffLengthX = this.boxSize.getXLength()/(THREE * (double) tmpCellNumberThirdX);
            // Multiple of 2 in y direction:
            int tmpCellNumberHalfY = (int) (this.boxSize.getYLength()/(TWO * this.cutOffLength));
            this.cellNumberY = INT_TWO * tmpCellNumberHalfY;
            if (this.cellNumberY%INT_TWO != 0) {
                throw new ArithmeticException("CellBox.initialise: this.cellNumberY can not be divided by 2 without rest");
            }
            this.cutOffLengthY = this.boxSize.getYLength()/(TWO * (double) tmpCellNumberHalfY);
            this.cellNumberXtimesY = this.cellNumberX * this.cellNumberY;
            // Multiple of 3 in x direction:
            int tmpCellNumberThirdZ = (int) (this.boxSize.getZLength()/(THREE * this.cutOffLength));
            this.cellNumberZ = INT_THREE * tmpCellNumberThirdZ;
            if (this.cellNumberZ%INT_THREE != 0) {
                throw new ArithmeticException("CellBox.initialise: this.cellNumberZ can not be divided by 3 without rest");
            }
            this.cutOffLengthZ = this.boxSize.getZLength()/(THREE * (double) tmpCellNumberThirdZ);
            this.cellNumber = this.cellNumberX * this.cellNumberY * this.cellNumberZ;
            this.cellChunkSize = this.cellNumberX/INT_THREE * this.cellNumberY/INT_TWO * this.cellNumberZ/INT_THREE;
            if (this.cellNumberX < INT_THREE || this.cellNumberY < INT_TWO || this.cellNumberZ < INT_THREE) {
                // There are less than minimum number of cells along a 
                // direction: NO cell-based calculation!
                this.isCellBasedCalculation = false;
                this.cellNeigbours = null;
                this.parallelizationSafeCellChunks = null;
                // <editor-fold defaultstate="collapsed" desc="Method call logging">
                this.simulationLogger.appendMethodCall("CellBox.Constructor: FULL, IsCellBasedCalculation = false");
                // </editor-fold>
            } else {
                this.isCellBasedCalculation = true;
                if (aCellNeigbours == null) {
                    this.cellNeigbours = this.getCellNeighbours();
                } else {
                    if (aCellNeigbours.length != this.cellNumber) {
                        throw new IllegalArgumentException("CellBox.Constructor: aCellNeigbours is illegal.");
                    }
                    this.cellNeigbours = aCellNeigbours;
                }
                if (aParallelizationSafeCellChunks == null) {
                    this.parallelizationSafeCellChunks = this.getParallelisationSafeCellChunks();
                } else {
                    this.parallelizationSafeCellChunks = aParallelizationSafeCellChunks;
                }
                // <editor-fold defaultstate="collapsed" desc="Method call logging">
                this.simulationLogger.appendMethodCall("CellBox.Constructor: FULL, IsCellBasedCalculation = true");
                // </editor-fold>
            }
        } catch (Exception anException) {
            // <editor-fold defaultstate="collapsed" desc="Exception logging">
            this.simulationLogger.appendException("CellBox.Constructor", Utils.getStacktrace(anException));
            // </editor-fold>
            throw anException;
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Protected methods">
    /**
     * Returns index of cell in array this.cells
     * NOTE: No checks are performed.
     * 
     * @param anIndexX Index in x-direction
     * @param anIndexY Index in x-direction
     * @param anIndexZ Index in x-direction
     * @return Index of cell in array this.cells
     */
    protected int getCellIndex(int anIndexX, int anIndexY, int anIndexZ) {
        // <editor-fold defaultstate="collapsed" desc="Explanation">
        // max(anIndexX) = this.cellNumberX - 1
        // max(anIndexY) = this.cellNumberY - 1
        // max(anIndexZ) = this.cellNumberZ - 1
        // 
        // anIndexX = 0; anIndexY = 0; anIndexZ = 0:
        // anCellIndex = 0
        // 
        // anIndexX = cellNumberX - 1; anIndexY = cellNumberY - 1; anIndexZ = cellNumberZ - 1:
        // cellIndex = cellNumberX - 1 + (cellNumberY - 1) *  cellNumberX + (cellNumberZ - 1) * cellNumberX * cellNumberY
        //           = cellNumberX - 1 + cellNumberY * cellNumberX - cellNumberX + cellNumberZ * cellNumberX * cellNumberY - cellNumberX * cellNumberY
        //           = - 1 + cellNumberY * cellNumberX + cellNumberZ * cellNumberX * cellNumberY - cellNumberX * cellNumberY
        //           = - 1 + cellNumberZ * cellNumberX * cellNumberY
        //           = cellNumberZ * cellNumberX * cellNumberY - 1
        // 
        // Correct:
        // return anIndexX + anIndexY * this.cellNumberX + anIndexZ * this.cellNumberX * this.cellNumberY;
        // but faster:
        // </editor-fold>
        return anIndexX + anIndexY * this.cellNumberX + anIndexZ * this.cellNumberXtimesY;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Returns cell neighbours
     * 
     * @return Cell neighbours
     */
    private int[][] getCellNeighbours() {
        // <editor-fold defaultstate="collapsed" desc="Explanation">
        // max(tmpIndexX) = this.cellNumberX - 1
        // max(tmpIndexY) = this.cellNumberY - 1
        // max(tmpIndexZ) = this.cellNumberZ - 1
        //
        // offset = cellNumberY * cellNumberX
        // 
        // -------------------------------------------------------------------
        // Lowest x,y plane:
        // -------------------------------------------------------------------
        // (cellNumberY - 1) * cellNumberX  ...  (cellNumberY - 1) * cellNumberX + cellNumberX - 1 = cellNumberY * cellNumberX - 1 
        //                                     = offset - 1
        // 
        // .                                     .
        // .                                     .
        // .                                     .
        // 
        // 2 * cellNumberX                  ...  3 * cellNumberX - 1
        // cellNumberX                      ...  2 * cellNumberX - 1
        // 0                                ...  cellNumberX - 1  
        // 
        // -------------------------------------------------------------------
        // Next upper x,y plane:
        // -------------------------------------------------------------------
        // offset + (cellNumberY - 1) * cellNumberX  ...  offset + (cellNumberY - 1) * cellNumberX + cellNumberX - 1 = offset + cellNumberY * cellNumberX - 1 
        //                                              = 2 * offset - 1
        // 
        // .                                              .
        // .                                              .
        // .                                              .
        // 
        // offset + 2 * cellNumberX                  ...  offset + 3 * cellNumberX - 1
        // offset + cellNumberX                      ...  offset + 2 * cellNumberX - 1
        // offset                                    ...  cellNumberX - 1  
        //
        // -------------------------------------------------------------------
        // Highest x,y plane:
        // -------------------------------------------------------------------
        // (cellNumberZ - 1) * offset + (cellNumberY - 1) * cellNumberX  ...  (cellNumberZ - 1) * offset + (cellNumberY - 1) * cellNumberX + cellNumberX - 1
        //                                                                  = (cellNumberZ - 1) * cellNumberY * cellNumberX + (cellNumberY - 1) * cellNumberX + cellNumberX - 1
        //                                                                  = cellNumberZ * cellNumberY * cellNumberX - cellNumberY * cellNumberX + cellNumberY * cellNumberX - cellNumberX + cellNumberX - 1
        //                                                                  = cellNumberZ * offset - 1
        //                                                                  = cellNumberZ * cellNumberY * cellNumberX - 1
        //
        // .                                                                  .
        // .                                                                  .
        // .                                                                  .
        // 
        // (cellNumberZ - 1) * offset + 2 * cellNumberX                  ...  (cellNumberZ - 1) * offset + 3 * cellNumberX - 1
        // (cellNumberZ - 1) * offset + cellNumberX                      ...  (cellNumberZ - 1) * offset + 2 * cellNumberX - 1
        // (cellNumberZ - 1) * offset                                    ...  (cellNumberZ - 1) * offset + cellNumberX - 1  
        //
        // </editor-fold>
        int[][] tmpCellNeigbours = new int[this.cellNumber][];
        for (int tmpIndexX = 0; tmpIndexX < this.cellNumberX; tmpIndexX++) {
            for (int tmpIndexY = 0; tmpIndexY < this.cellNumberY; tmpIndexY++) {
                for (int tmpIndexZ = 0; tmpIndexZ < this.cellNumberZ; tmpIndexZ++) {
                    // There are 13 neighbour cells (without double counting of later pair interactions)
                    int tmpCurrentCellIndex = this.getCellIndex(tmpIndexX, tmpIndexY, tmpIndexZ);
                    tmpCellNeigbours[tmpCurrentCellIndex] = new int[13];
                    Arrays.fill(tmpCellNeigbours[tmpCurrentCellIndex], -1);
                    int tmpNeighbourIndex = 0;
                    // <editor-fold defaultstate="collapsed" desc="1. Neighbour cell: (x + 1, y, z)">
                    int tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    int tmpCorrectedIndexY = this.correctIndexY(tmpIndexY);
                    int tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="2. Neighbour cell: (x + 1, y + 1, z)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="3. Neighbour cell: (x, y + 1, z)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="4. Neighbour cell: (x - 1, y + 1, z)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX - 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="5. Neighbour cell: (x + 1, y, z - 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ - 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="6. Neighbour cell: (x + 1, y + 1, z - 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ - 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="7. Neighbour cell: (x, y + 1, z - 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ - 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="8. Neighbour cell: (x - 1, y + 1, z - 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX - 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ - 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="9. Neighbour cell: (x + 1, y, z + 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ + 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="10. Neighbour cell: (x + 1, y + 1, z + 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX + 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ + 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="11. Neighbour cell: (x, y + 1, z + 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ + 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="12. Neighbour cell: (x - 1, y + 1, z + 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX - 1);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY + 1);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ + 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                    // <editor-fold defaultstate="collapsed" desc="13. Neighbour cell: (x, y, z + 1)">
                    tmpCorrectedIndexX = this.correctIndexX(tmpIndexX);
                    tmpCorrectedIndexY = this.correctIndexY(tmpIndexY);
                    tmpCorrectedIndexZ = this.correctIndexZ(tmpIndexZ + 1);
                    if (tmpCorrectedIndexX >= 0 && tmpCorrectedIndexY >= 0 && tmpCorrectedIndexZ >= 0) {
                        int tmpNeighbourCellIndex = this.getCellIndex(tmpCorrectedIndexX, tmpCorrectedIndexY, tmpCorrectedIndexZ);
                        tmpCellNeigbours[tmpCurrentCellIndex][tmpNeighbourIndex++] = tmpNeighbourCellIndex;
                    }
                    // </editor-fold>
                }                
            }
        }
        return tmpCellNeigbours;
    }

    /**
     * Returns parallelisation-safe cell chunks
     * 
     * @return Parallelisation-safe cell chunks
     */
    private int[][] getParallelisationSafeCellChunks() {
        // Dimension 18: 3 (in x direction) * 2 (in y direction) * 3 (in z direction) = 18
        int[][] tmpParallelisationSafeCellChunks = new int[18][this.cellChunkSize];
        int tmpChunkIndex = 0;
        for (int tmpOffsetZ = 0; tmpOffsetZ < 3; tmpOffsetZ++) {
            for (int tmpOffsetY = 0; tmpOffsetY < 2; tmpOffsetY++) {
                for (int tmpOffsetX = 0; tmpOffsetX < 3; tmpOffsetX++) {
                    int tmpCellIndex = 0;
                    for (int tmpIndexZ = 0; tmpIndexZ < this.cellNumberZ - 1; tmpIndexZ += 3) {
                        for (int tmpIndexY = 0; tmpIndexY < this.cellNumberY - 1; tmpIndexY += 2) {
                            for (int tmpIndexX= 0; tmpIndexX < this.cellNumberX - 1; tmpIndexX += 3) {
                                tmpParallelisationSafeCellChunks[tmpChunkIndex][tmpCellIndex++] = 
                                    this.getCellIndex(tmpOffsetX + tmpIndexX, tmpOffsetY + tmpIndexY, tmpOffsetZ + tmpIndexZ);
                            }
                        }
                    }
                    tmpChunkIndex++;
                }
            }
        }
        return tmpParallelisationSafeCellChunks;
    }
    
    /**
     * Corrects index in x-direction
     * 
     * @param anIndexX Index in x-direction
     * @return Corrected index or -1 if index is illegal
     */
    private int correctIndexX(int anIndexX) {
        if (anIndexX < 0) {
            if (this.periodicBoundaries.isAlongX()) {
                return this.cellNumberX - 1;
            } else {
                return -1;
            }
        } else if (anIndexX == this.cellNumberX) {
            if (this.periodicBoundaries.isAlongX()) {
                return 0;
            } else {
                return -1;
            }
        } else {
            return anIndexX;
        }
    }
    
    /**
     * Corrects index in y-direction
     * 
     * @param anIndexY Index in y-direction
     * @return Corrected index or -1 if index is illegal
     */
    private int correctIndexY(int anIndexY) {
        if (anIndexY < 0) {
            if (this.periodicBoundaries.isAlongY()) {
                return this.cellNumberY - 1;
            } else {
                return -1;
            }
        } else if (anIndexY == this.cellNumberY) {
            if (this.periodicBoundaries.isAlongY()) {
                return 0;
            } else {
                return -1;
            }
        } else {
            return anIndexY;
        }
    }

    /**
     * Corrects index in z-direction
     * 
     * @param anIndexZ Index in z-direction
     * @return Corrected index or -1 if index is illegal
     */
    private int correctIndexZ(int anIndexZ) {
        if (anIndexZ < 0) {
            if (this.periodicBoundaries.isAlongZ()) {
                return this.cellNumberZ - 1;
            } else {
                return -1;
            }
        } else if (anIndexZ == this.cellNumberZ) {
            if (this.periodicBoundaries.isAlongZ()) {
                return 0;
            } else {
                return -1;
            }
        } else {
            return anIndexZ;
        }
    }
    // </editor-fold>

}
