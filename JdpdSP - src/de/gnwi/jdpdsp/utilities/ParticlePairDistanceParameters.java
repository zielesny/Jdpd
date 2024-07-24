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

import java.util.LinkedList;

/**
 * Particle pair distance parameters
 * 
 * @author Achim Zielesny
 */
public class ParticlePairDistanceParameters {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Capacity (initial and growth increment)
     */
    private final int capacity;
    
    /**
     * List for particle i index arrays
     */
    private final LinkedList<int[]> particleIndex_i_ArrayList;
    
    /**
     * List for particle j index arrays
     */
    private final LinkedList<int[]> particleIndex_j_ArrayList;

    /**
     * List for Rij_x arrays
     */
    private final LinkedList<float[]> arrayList_Rij_x;

    /**
     * List for Rij_y arrays
     */
    private final LinkedList<float[]> arrayList_Rij_y;

    /**
     * List for Rij_z arrays
     */
    private final LinkedList<float[]> arrayList_Rij_z;

    /**
     * List for Rij_Square arrays
     */
    private final LinkedList<float[]> arrayList_Rij_Square;

    /**
     * List for Rij arrays
     */
    private final LinkedList<float[]> arrayList_Rij;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Particle i indices
     */
    private int[] particleIndices_i;
    
    /**
     * Particle j indices
     */
    private int[] particleIndices_j;
    
    /**
     * Rij_x array
     */
    private float[] array_Rij_x;
    
    /**
     * Rij_y array
     */
    private float[] array_Rij_y;
    
    /**
     * Rij_z array
     */
    private float[] array_Rij_z;
    
    /**
     * Rij_Square array
     */
    private float[] array_Rij_Square;
    
    /**
     * Rij array
     */
    private float[] array_Rij;
    
    /**
     * Index
     */
    private int index;
    
    /**
     * Size (exclusive end index)
     */
    private int size;
    
    /**
     * Maximum index
     */
    private int maxIndex;
    
    /**
     * True: Arrays are consolidated, false: Otherwise
     */
    private boolean isConsolidated;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aCapacity Capacity (initial and growth increment)
     * @throws IllegalArgumentException Thrown if an argument is illegal
     */
    public ParticlePairDistanceParameters(int aCapacity) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aCapacity < 1) {
            throw new IllegalArgumentException("ParticlePairInfoAdder.Constructor: aCapacity < 1.");
        }
        // </editor-fold>
        this.capacity = aCapacity;
        this.particleIndex_i_ArrayList = new LinkedList<>();
        this.particleIndex_j_ArrayList = new LinkedList<>();
        this.arrayList_Rij_x = new LinkedList<>();
        this.arrayList_Rij_y = new LinkedList<>();
        this.arrayList_Rij_z = new LinkedList<>();
        this.arrayList_Rij_Square = new LinkedList<>();
        this.arrayList_Rij = new LinkedList<>();

        this.particleIndices_i = new int[this.capacity];
        this.particleIndices_j = new int[this.capacity];
        this.array_Rij_x = new float[this.capacity];
        this.array_Rij_y = new float[this.capacity];
        this.array_Rij_z = new float[this.capacity];
        this.array_Rij_Square = new float[this.capacity];
        this.array_Rij = new float[this.capacity];
        
        this.maxIndex = this.capacity - 1;

        this.reset();
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Adds particle pair parameters
     * 
     * @param aParticleIndex_i Index of particle i
     * @param aParticleIndex_j Index of particle j
     * @param aRij_x Rij_x
     * @param aRij_y Rij_y
     * @param aRij_z Rij_z
     * @param aRij_Square Rij_square
     * @param aRij Rij
     */
    public void add(
        int aParticleIndex_i, 
        int aParticleIndex_j,
        float aRij_x,
        float aRij_y,
        float aRij_z,
        float aRij_Square,
        float aRij
    ) {
        if (this.index > this.maxIndex) {
            // IMPORTANT: Consolidation status is false
            this.isConsolidated = false;

            this.particleIndex_i_ArrayList.add(this.particleIndices_i);
            this.particleIndex_j_ArrayList.add(this.particleIndices_j);
            this.arrayList_Rij_x.add(this.array_Rij_x);
            this.arrayList_Rij_y.add(this.array_Rij_y);
            this.arrayList_Rij_z.add(this.array_Rij_z);
            this.arrayList_Rij_Square.add(this.array_Rij_Square);
            this.arrayList_Rij.add(this.array_Rij);

            this.particleIndices_i = new int[this.capacity];
            this.particleIndices_j = new int[this.capacity];
            this.array_Rij_x = new float[this.capacity];
            this.array_Rij_y = new float[this.capacity];
            this.array_Rij_z = new float[this.capacity];
            this.array_Rij_Square = new float[this.capacity];
            this.array_Rij = new float[this.capacity];

            this.maxIndex = this.capacity - 1;
            this.index = 0;
        }
        this.particleIndices_i[this.index] = aParticleIndex_i;
        this.particleIndices_j[this.index] = aParticleIndex_j;
        this.array_Rij_x[this.index] = aRij_x;
        this.array_Rij_y[this.index] = aRij_y;
        this.array_Rij_z[this.index] = aRij_z;
        this.array_Rij_Square[this.index] = aRij_Square;
        this.array_Rij[this.index] = aRij;
        
        this.index++;
        this.size++;
    }
    
    /**
     * Resets particle pair parameters
     */
    public final void reset() {
        this.particleIndex_i_ArrayList.clear();
        this.particleIndex_j_ArrayList.clear();
        this.arrayList_Rij_x.clear();
        this.arrayList_Rij_y.clear();
        this.arrayList_Rij_z.clear();
        this.arrayList_Rij_Square.clear();
        this.arrayList_Rij.clear();
        
        this.index = 0;
        this.size = 0;
        
        // IMPORTANT: After reset the consolidation status is true 
        this.isConsolidated = true;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Size (exclusive end index)
     * 
     * @return Size (exclusive end index)
     */
    public int getSize() {
        return this.size;
    }

    /**
     * Particle i indices
     * 
     * @return Particle i indices
     * @throws IllegalStateException Thrown if arrays are not consolidated
     */
    public int[] getParticleIndices_i() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.particleIndices_i;
    }

    /**
     * Particle j indices
     * 
     * @return Particle j indices
     */
    public int[] getParticleIndices_j() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.particleIndices_j;
    }

    /**
     * Rij_x array
     * 
     * @return Rij_x array
     */
    public float[] getArray_Rij_x() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.array_Rij_x;
    }

    /**
     * Rij_y array
     * 
     * @return Rij_y array
     */
    public float[] getArray_Rij_y() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.array_Rij_y;
    }

    /**
     * Rij_z array
     * 
     * @return Rij_z array
     */
    public float[] getArray_Rij_z() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.array_Rij_z;
    }

    /**
     * Rij_Square array
     * 
     * @return Rij_Square array
     */
    public float[] getArray_Rij_Square() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.array_Rij_Square;
    }

    /**
     * Rij array
     * 
     * @return Rij array
     */
    public float[] getArray_Rij() {
        if (!this.isConsolidated) {
            this.consolidate();
        }
        return this.array_Rij;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    
    /**
     * Consolidates arrays
     * (No checks are performed)
     */
    private void consolidate() {
        if (!this.isConsolidated) {
            int tmpAccumulatedArraySize = 0;
            for (int[] tmpParticleIndexArray : this.particleIndex_i_ArrayList) {
                tmpAccumulatedArraySize += tmpParticleIndexArray.length;
            }
            int tmpNewArraySize = tmpAccumulatedArraySize + this.capacity;
            this.maxIndex = tmpNewArraySize - 1;

            this.particleIndices_i = 
                this.getConsolidatedIntegerArray(tmpNewArraySize,
                    this.particleIndex_i_ArrayList,
                    this.particleIndices_i
                );
            this.particleIndices_j = 
                this.getConsolidatedIntegerArray(tmpNewArraySize,
                    this.particleIndex_j_ArrayList,
                    this.particleIndices_j
                );
            this.array_Rij_x = 
                this.getConsolidatedFloatArray(
                    tmpNewArraySize,
                    this.arrayList_Rij_x,
                    this.array_Rij_x
                );
            this.array_Rij_y = 
                this.getConsolidatedFloatArray(
                    tmpNewArraySize,
                    this.arrayList_Rij_y,
                    this.array_Rij_y
                );
            this.array_Rij_z = 
                this.getConsolidatedFloatArray(
                    tmpNewArraySize,
                    this.arrayList_Rij_z,
                    this.array_Rij_z
                );
            this.array_Rij_Square = 
                this.getConsolidatedFloatArray(
                    tmpNewArraySize,
                    this.arrayList_Rij_Square,
                    this.array_Rij_Square
                );
            this.array_Rij = 
                this.getConsolidatedFloatArray(tmpNewArraySize,
                    this.arrayList_Rij,
                    this.array_Rij
                );
            
            this.particleIndex_i_ArrayList.clear();
            this.particleIndex_j_ArrayList.clear();
            this.arrayList_Rij_x.clear();
            this.arrayList_Rij_y.clear();
            this.arrayList_Rij_z.clear();
            this.arrayList_Rij_Square.clear();
            this.arrayList_Rij.clear();
            
            this.isConsolidated = true;
        }
    }
    
    /**
     * Consolidates integer arrays
     * (No checks are performed)
     * 
     * @param tmpNewArraySize New array size
     * @param anArrayList Array list with arrays to consolidate
     * @param aLastArray Last array to be consolidated
     * @return Consolidated integer array
     */
    private int[] getConsolidatedIntegerArray(
        int tmpNewArraySize,
        LinkedList<int[]> anArrayList,
        int[] aLastArray
    ) {
        int[] tmpNewArray = new int[tmpNewArraySize];
        int tmpIndex = 0;
        for (int[] tmpSingleArray : anArrayList) {
            System.arraycopy(
                tmpSingleArray, 
                0, 
                tmpNewArray, 
                tmpIndex, 
                tmpSingleArray.length
            );
            tmpIndex += tmpSingleArray.length;
        }
        System.arraycopy(
            aLastArray, 
            0, 
            tmpNewArray, 
            tmpIndex, 
            aLastArray.length
        );
        return tmpNewArray;
    }

    /**
     * Consolidates float arrays
     * (No checks are performed)
     * 
     * @param tmpNewArraySize New array size
     * @param anArrayList Array list with arrays to consolidate
     * @param aLastArray Last array to be consolidated
     * @return Consolidated float array
     */
    private float[] getConsolidatedFloatArray(
        int tmpNewArraySize,
        LinkedList<float[]> anArrayList,
        float[] aLastArray
    ) {
        float[] tmpNewArray = new float[tmpNewArraySize];
        int tmpIndex = 0;
        for (float[] tmpSingleArray : anArrayList) {
            System.arraycopy(
                tmpSingleArray, 
                0, 
                tmpNewArray, 
                tmpIndex, 
                tmpSingleArray.length
            );
            tmpIndex += tmpSingleArray.length;
        }
        System.arraycopy(
            aLastArray, 
            0, 
            tmpNewArray, 
            tmpIndex, 
            aLastArray.length
        );
        return tmpNewArray;
    }
    // </editor-fold>
    
}
