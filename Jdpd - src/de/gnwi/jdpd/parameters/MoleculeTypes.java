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

import java.util.HashMap;

/**
 * Molecule types
 * 
 * @author Achim Zielesny
 */
public class MoleculeTypes {
    
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Number of molecule types
     */
    private final int moleculeTypeNumber;
    
    /**
     * Map of molecule name to molecule index. This index refers to molecule related arrays.
     */
    private final HashMap<String, Integer> moleculeNameToIndexMap;
    
    /**
     * Molecule names
     */
    private final String[] moleculeNames;
    
    /**
     * Molar masses of molecules
     */
    private final double[] moleculeMolarMasses;
    
    /**
     * Total number of all particles of molecules
     */
    private final int[] totalMoleculeParticleNumbers;

    /**
     * Number of particles of a single molecule
     */
    private final int[] singleMoleculeParticleNumbers;
    
    /**
     * First indices in particle arrays
     */
    private final int[] firstIndices;
    
    /**
     * Last indices in particle arrays
     */
    private final int[] lastIndices;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)

     * @param aMoleculeNameToIndexMap Map of molecule name to molecule index. This index refers to molecule related arrays.
     * @param aMoleculeNames Molecule names
     * @param aMoleculeMolarMasses Molar masses of molecules
     * @param aTotalMoleculeParticleNumbers Total number of all particles of molecules
     * @param aSingleMoleculeParticleNumbers Number of particles of a single molecule
     * @param aFirstIndices First indices in particle arrays
     * @param aLastIndices Last indices in particle arrays
     */
    public MoleculeTypes(
        HashMap<String, Integer> aMoleculeNameToIndexMap,
        String[] aMoleculeNames,
        double[] aMoleculeMolarMasses,
        int[] aTotalMoleculeParticleNumbers,
        int[] aSingleMoleculeParticleNumbers,
        int[] aFirstIndices,
        int[] aLastIndices) {
        this.moleculeTypeNumber = aMoleculeNameToIndexMap.size();
        this.moleculeNameToIndexMap = aMoleculeNameToIndexMap;
        this.moleculeNames = aMoleculeNames;
        this.moleculeMolarMasses = aMoleculeMolarMasses;
        this.totalMoleculeParticleNumbers = aTotalMoleculeParticleNumbers;
        this.singleMoleculeParticleNumbers = aSingleMoleculeParticleNumbers;
        this.firstIndices = aFirstIndices;
        this.lastIndices = aLastIndices;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Number of molecule types
     * 
     * @return Number of molecule types
     */
    public int getMoleculeTypeNumber() {
        return this.moleculeTypeNumber;
    }
    
    /**
     * Index of molecule name
     * 
     * @param aMoleculeName MoleculeName
     * @return Index of molecule name
     */
    public int getIndex(String aMoleculeName) {
        return this.moleculeNameToIndexMap.get(aMoleculeName);
    }
    
    /**
     * Molecule names
     * 
     * @return Molecule names
     */
    public String[] getMoleculeNames() {
        return this.moleculeNames;
    }
    
    /**
     * Molecule name
     * 
     * @param anIndex Molecule index
     * @return Molecule name
     */
    public String getMoleculeName(int anIndex) {
        return this.moleculeNames[anIndex];
    }
    
    /**
     * Total number of all particles of molecules
     * 
     * @param aMoleculeName Molecule name
     * @return Total number of all particles of molecules
     */
    public int getTotalMoleculeParticleNumber(String aMoleculeName) {
        return this.totalMoleculeParticleNumbers[this.moleculeNameToIndexMap.get(aMoleculeName)];
    }
    
    /**
     * Total number of all particles of molecules
     * 
     * @param anIndex Molecule index
     * @return Total number of all particles of molecules
     */
    public int getTotalMoleculeParticleNumber(int anIndex) {
        return this.totalMoleculeParticleNumbers[anIndex];
    }

    /**
     * Number of particles of a single molecule
     * 
     * @param aMoleculeName Molecule name
     * @return Number of particles of a single molecule
     */
    public int getSingleMoleculeParticleNumber(String aMoleculeName) {
        return this.singleMoleculeParticleNumbers[this.moleculeNameToIndexMap.get(aMoleculeName)];
    }

    /**
     * Number of particles of a single molecule
     * 
     * @param anIndex Molecule index
     * @return Number of particles of a single molecule
     */
    public int getSingleMoleculeParticleNumber(int anIndex) {
        return this.singleMoleculeParticleNumbers[anIndex];
    }
    
    /**
     * First index in particle arrays
     * 
     * @param aMoleculeName Molecule name
     * @return First index in particle arrays
     */
    public int getFirstIndex(String aMoleculeName) {
        return this.firstIndices[this.moleculeNameToIndexMap.get(aMoleculeName)];
    }
    
    /**
     * First index in particle arrays
     * 
     * @param anIndex Molecule index
     * @return First index in particle arrays
     */
    public int getFirstIndex(int anIndex) {
        return this.firstIndices[anIndex];
    }
    
    /**
     * Last index in particle arrays
     * 
     * @param aMoleculeName Molecule name
     * @return Last index in particle arrays
     */
    public int getLastIndex(String aMoleculeName) {
        return this.lastIndices[this.moleculeNameToIndexMap.get(aMoleculeName)];
    }
    
    /**
     * Last index in particle arrays
     * 
     * @param anIndex Molecule index
     * @return Last index in particle arrays
     */
    public int getLastIndex(int anIndex) {
        return this.lastIndices[anIndex];
    }

    /**
     * Molar mass of molecule
     * 
     * @param aMoleculeName Molecule name
     * @return Molar mass of molecule
     */
    public double getMoleculeMolarMass(String aMoleculeName) {
        return this.moleculeMolarMasses[this.moleculeNameToIndexMap.get(aMoleculeName)];
    }
    
    /**
     * Molar mass of molecule
     * 
     * @param anIndex Molecule index
     * @return Molar mass of molecule
     */
    public double getMoleculeMolarMass(int anIndex) {
        return this.moleculeMolarMasses[anIndex];
    }
    // </editor-fold>
    
}
