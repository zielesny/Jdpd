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
package de.gnwi.jdpd.nearestNeighbor;

import de.gnwi.jdpd.utilities.Utils;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Set;

/**
 * Nearest-neighbor manager
 * 
 * @author Achim Zielesny
 */
public class NearestNeighborManager {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * (Non-defined) Index value that indicates no nearest neighbor
     */
    private static final int NO_NEAREST_NEIGHBOR_INDEX_VALUE = -1;
    
    /**
     * String that indicates no nearest neighbor
     */
    private static final String NO_NEAREST_NEIGHBOR_STRING = "None";
    
    /**
     * Molecule name separator
     */
    private static final String MOLECULE_NAME_SEPARATOR = ", ";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * List with all base particle indices in particle arrays for which the 
     * nearest-neighbor determination is to be performed
     */
    private final LinkedList<Integer> baseParticleIndexList;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private class variables">
    /**
     * Array with all base particle indices in particle arrays for which the 
     * nearest-neighbor determination is to be performed
     */
    private int[] baseParticleIndices;

    /**
     * Base molecule-particle to nearest-neighbor molecule-particle frequency map
     */
    private HashMap<String, HashMap<String, Integer>> baseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap;

    /**
     * Base molecule-particle to nearest-neighbor particle frequency map
     */
    private HashMap<String, HashMap<String, Integer>> baseMoleculeParticleToNearestNeighborParticleFrequencyMap;

    /**
     * Base molecule-particle to nearest-neighbor molecule frequency map
     */
    private HashMap<String, HashMap<String, Integer>> baseMoleculeParticleToNearestNeighborMoleculeFrequencyMap;

    /**
     * Base molecule to nearest-neighbor molecule frequency map
     */
    private HashMap<String, HashMap<String, Integer>> baseMoleculeToNearestNeighborMoleculeFrequencyMap;

    /**
     * Base molecule to nearest-neighbor molecule-tuple frequency map
     */
    private HashMap<String, HashMap<String, Integer>> baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap;
    
    /**
     * Nearest-neighbor molecule name map
     */
    private HashMap<String, String> nearestNeighborMoleculeNameMap;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * NOTE: NO checks are performed
     */
    public NearestNeighborManager() {
        this.baseParticleIndexList = new LinkedList<>();
        this.baseParticleIndices = null;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Adds NearestNeighborBaseParticleDescription instance
     * 
     * @param aNearestNeighborBaseParticleDescription NearestNeighborBaseParticleDescription instance
     */
    public void addNearestNeighborParticleDescription(NearestNeighborBaseParticleDescription aNearestNeighborBaseParticleDescription) {
        this.baseParticleIndexList.addAll(aNearestNeighborBaseParticleDescription.getBaseParticleIndexList());
    }
    
    /**
     * Consolidates this.baseParticleIndexList to array this.baseParticleIndices
     */
    public void consolidate() {
        this.baseParticleIndices = new int[this.baseParticleIndexList.size()];
        int tmpIndex = 0;
        for (Integer tmpBaseParticleIndex : this.baseParticleIndexList) {
            this.baseParticleIndices[tmpIndex++] = tmpBaseParticleIndex;
        }
        // IMPORTANT: SORT this.baseParticleIndices to correctly evaluate 
        //            particles that belong to the same molecule (which must be
        //            consecutively one after another)
        Arrays.sort(this.baseParticleIndices);
        this.baseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap = new HashMap<>(this.baseParticleIndices.length);
        this.baseMoleculeParticleToNearestNeighborParticleFrequencyMap = new HashMap<>(this.baseParticleIndices.length);
        this.baseMoleculeParticleToNearestNeighborMoleculeFrequencyMap = new HashMap<>(this.baseParticleIndices.length);
        // Hash maps will be very small for most practical purposes so default hash map size will be sufficient
        this.baseMoleculeToNearestNeighborMoleculeFrequencyMap = new HashMap<>();
        this.baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap = new HashMap<>();
        // Hash map will be very small for most practical purposes so default hash map size will be sufficient
        this.nearestNeighborMoleculeNameMap = new HashMap<>();
    }
    
    /**
     * Re-initializes specified arrays.
     * NOTE: NO checks are performed
     * 
     * @param aNearestNeighborDistances Distances to nearest (non-bonded) 
     * neighbor particles (of another molecule)
     * @param aNearestNeighborParticleIndices Indices of nearest (non-bonded) 
     * neighbor particles (of another molecule)
     */
    public void reInitialize(
        double[] aNearestNeighborDistances,
        int[] aNearestNeighborParticleIndices
    ) {
        for (int tmpBaseParticleIndex : this.baseParticleIndices) {
            aNearestNeighborDistances[tmpBaseParticleIndex] = Double.MAX_VALUE;
            aNearestNeighborParticleIndices[tmpBaseParticleIndex] = NearestNeighborManager.NO_NEAREST_NEIGHBOR_INDEX_VALUE;
        }
    }

    /**
     * Evaluates nearest-neighbors
     * NOTE: NO checks are performed.
     * 
     * @param aNearestNeighborParticleIndices Nearest-neighbor particle indices in particle arrays
     * @param aParticleTokens Particle tokens in particle arrays
     * @param aMoleculeNamesOfParticles Molecule names of particles in particle arrays 
     * @param aMoleculeIndices Indeces of molecules
     */
    public void evaluateNearestNeighbors(
        int[] aNearestNeighborParticleIndices,
        String[] aParticleTokens,
        String[] aMoleculeNamesOfParticles,
        int[] aMoleculeIndices
    ) {
        this.baseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap.clear();
        this.baseMoleculeParticleToNearestNeighborParticleFrequencyMap.clear();
        this.baseMoleculeParticleToNearestNeighborMoleculeFrequencyMap.clear();
        this.baseMoleculeToNearestNeighborMoleculeFrequencyMap.clear();
        this.baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap.clear();
        this.nearestNeighborMoleculeNameMap.clear();
        int tmpOldBaseMoleculeIndex = -1;
        String tmpOldBaseMoleculeName = null;
        for (int i = 0; i < this.baseParticleIndices.length; i++) {
            int tmpBaseParticleIndex = this.baseParticleIndices[i];
            // Base particle
            String tmpBaseParticleToken = aParticleTokens[tmpBaseParticleIndex];
            String tmpBaseMoleculeName = aMoleculeNamesOfParticles[tmpBaseParticleIndex];
            String tmpBaseMoleculeParticle = Utils.getMoleculeParticle(tmpBaseParticleToken, tmpBaseMoleculeName);
            int tmpBaseMoleculeIndex = aMoleculeIndices[tmpBaseParticleIndex];
            // Nearest-neighbor particle
            String tmpNearestNeighborParticleToken = NearestNeighborManager.NO_NEAREST_NEIGHBOR_STRING;
            String tmpNearestNeighborMoleculeName = NearestNeighborManager.NO_NEAREST_NEIGHBOR_STRING;
            String tmpNearestNeighborMoleculeParticle = NearestNeighborManager.NO_NEAREST_NEIGHBOR_STRING;
            int tmpNearestNeighborMoleculeIndex = NearestNeighborManager.NO_NEAREST_NEIGHBOR_INDEX_VALUE;
            int tmpNearestNeighborParticleIndex = aNearestNeighborParticleIndices[tmpBaseParticleIndex];
            if (tmpNearestNeighborParticleIndex != NearestNeighborManager.NO_NEAREST_NEIGHBOR_INDEX_VALUE) {
                tmpNearestNeighborParticleToken = aParticleTokens[tmpNearestNeighborParticleIndex];
                tmpNearestNeighborMoleculeName = aMoleculeNamesOfParticles[tmpNearestNeighborParticleIndex];
                tmpNearestNeighborMoleculeParticle = Utils.getMoleculeParticle(tmpNearestNeighborParticleToken, tmpNearestNeighborMoleculeName);
                tmpNearestNeighborMoleculeIndex = aMoleculeIndices[tmpNearestNeighborParticleIndex];
            }
            // Set molecule-particle nearest-neighbor hash maps
            this.setBaseToNearestNeighborFrequencyMap(
                tmpBaseMoleculeParticle,
                this.baseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap,
                tmpNearestNeighborMoleculeParticle
            );
            this.setBaseToNearestNeighborFrequencyMap(
                tmpBaseMoleculeParticle,
                this.baseMoleculeParticleToNearestNeighborParticleFrequencyMap,
                tmpNearestNeighborParticleToken
            );
            this.setBaseToNearestNeighborFrequencyMap(
                tmpBaseMoleculeParticle,
                this.baseMoleculeParticleToNearestNeighborMoleculeFrequencyMap,
                tmpNearestNeighborMoleculeName
            );
            // Set molecule nearest-neighbor hash map
            if (i == 0) {
                tmpOldBaseMoleculeIndex = tmpBaseMoleculeIndex;
                tmpOldBaseMoleculeName = tmpBaseMoleculeName;
                this.nearestNeighborMoleculeNameMap.put(tmpNearestNeighborMoleculeName, tmpNearestNeighborMoleculeName);
            } else {
                if (tmpOldBaseMoleculeIndex == tmpBaseMoleculeIndex) {
                    if (!this.nearestNeighborMoleculeNameMap.containsKey(tmpNearestNeighborMoleculeName)) {
                        this.nearestNeighborMoleculeNameMap.put(tmpNearestNeighborMoleculeName, tmpNearestNeighborMoleculeName);
                    }
                } else {
                    for (String tmpNearestNeighborMoleculeNameInMap : this.nearestNeighborMoleculeNameMap.keySet()) {
                        this.setBaseToNearestNeighborFrequencyMap(
                            tmpOldBaseMoleculeName,
                            this.baseMoleculeToNearestNeighborMoleculeFrequencyMap,
                            tmpNearestNeighborMoleculeNameInMap
                        );
                    }
                    String tmpNearestNeighborMoleculeNameTuple = this.getMoleculeNameTuple(this.nearestNeighborMoleculeNameMap.keySet());
                    this.setBaseToNearestNeighborFrequencyMap(
                        tmpOldBaseMoleculeName,
                        this.baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
                        tmpNearestNeighborMoleculeNameTuple
                    );
                    tmpOldBaseMoleculeIndex = tmpBaseMoleculeIndex;
                    tmpOldBaseMoleculeName = tmpBaseMoleculeName;
                    this.nearestNeighborMoleculeNameMap.clear();
                    this.nearestNeighborMoleculeNameMap.put(tmpNearestNeighborMoleculeName, tmpNearestNeighborMoleculeName);
                }
            }
            if (i == this.baseParticleIndices.length - 1) {
                for (String tmpNearestNeighborMoleculeNameInMap : this.nearestNeighborMoleculeNameMap.keySet()) {
                    this.setBaseToNearestNeighborFrequencyMap(
                        tmpOldBaseMoleculeName,
                        this.baseMoleculeToNearestNeighborMoleculeFrequencyMap,
                        tmpNearestNeighborMoleculeNameInMap
                    );
                }
                String tmpNearestNeighborMoleculeNameTuple = this.getMoleculeNameTuple(this.nearestNeighborMoleculeNameMap.keySet());
                this.setBaseToNearestNeighborFrequencyMap(
                    tmpOldBaseMoleculeName,
                    this.baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
                    tmpNearestNeighborMoleculeNameTuple
                );
            }
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Base molecule-particle to nearest-neighbor molecule-particle frequency map
     * 
     * @return Base molecule-particle to nearest-neighbor molecule-particle frequency map
     */
    public HashMap<String, HashMap<String, Integer>> getBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap() {
        return this.baseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap;
    }

    /**
     * Base molecule-particle to nearest-neighbor particle frequency map
     * 
     * @return Base molecule-particle to nearest-neighbor particle frequency map
     */
    public HashMap<String, HashMap<String, Integer>> getBaseMoleculeParticleToNearestNeighborParticleFrequencyMap() {
        return this.baseMoleculeParticleToNearestNeighborParticleFrequencyMap;
    }

    /**
     * Base molecule-particle to nearest-neighbor molecule frequency map
     * 
     * @return Base molecule-particle to nearest-neighbor molecule frequency map
     */
    public HashMap<String, HashMap<String, Integer>> getBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap() {
        return this.baseMoleculeParticleToNearestNeighborMoleculeFrequencyMap;
    }

    /**
     * Base molecule to nearest-neighbor molecule frequency map
     * 
     * @return Base molecule to nearest-neighbor molecule frequency map
     */
    public HashMap<String, HashMap<String, Integer>> getBaseMoleculeToNearestNeighborMoleculeFrequencyMap() {
        return this.baseMoleculeToNearestNeighborMoleculeFrequencyMap;
    }

    /**
     * Base molecule to nearest-neighbor molecule-tuple frequency map
     * 
     * @return Base molecule to nearest-neighbor molecule-tuple frequency map
     */
    public HashMap<String, HashMap<String, Integer>> getBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap() {
        return this.baseMoleculeToNearestNeighborMoleculeTupleFrequencyMap;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Sets aBaseToNearestNeighborFrequencyMap (see code)
     * NOTE: NO checks are performed
     * 
     * @param aBase Base
     * @param aBaseToNearestNeighborFrequencyMap Base to nearest-neighbor 
     * frequency map
     * @param aNearestNeighbor Nearest-neighbor
     */
    private void setBaseToNearestNeighborFrequencyMap(
        String aBase,
        HashMap<String, HashMap<String, Integer>> aBaseToNearestNeighborFrequencyMap,
        String aNearestNeighbor
    ) {
        if (aBaseToNearestNeighborFrequencyMap.containsKey(aBase)) {
            HashMap<String, Integer> tmpNearestNeighborToFrequencyMap = aBaseToNearestNeighborFrequencyMap.get(aBase);
            if (tmpNearestNeighborToFrequencyMap.containsKey(aNearestNeighbor)) {
                int tmpCurrentFrequency = tmpNearestNeighborToFrequencyMap.get(aNearestNeighbor);
                tmpNearestNeighborToFrequencyMap.replace(aNearestNeighbor, tmpCurrentFrequency + 1);
            } else {
                tmpNearestNeighborToFrequencyMap.put(aNearestNeighbor, 1);
            }
        } else {
            // Default hash map size is sufficient for most practical purposes
            HashMap<String, Integer> tmpNearestNeighborToFrequencyMap = new HashMap<>();
            tmpNearestNeighborToFrequencyMap.put(aNearestNeighbor, 1);
            aBaseToNearestNeighborFrequencyMap.put(aBase, tmpNearestNeighborToFrequencyMap);
        }
    }
    
    /**
     * Returns molecule name tuple
     * NOTE: NO checks are performed.
     * 
     * @param aMoleculeNameSet Molecule name set to generate tuple from
     * @return Molecule name tuple
     */
    private String getMoleculeNameTuple(Set<String> aMoleculeNameSet) {
        StringBuilder tmpTupleBuffer;
        switch (aMoleculeNameSet.size()) {
            case 1:
                return aMoleculeNameSet.iterator().next();
            case 2:
                Iterator<String> tmpMoleculeIterator = aMoleculeNameSet.iterator();
                String tmpMoleculeName1 = tmpMoleculeIterator.next();
                String tmpMoleculeName2 = tmpMoleculeIterator.next();
                tmpTupleBuffer = new StringBuilder(tmpMoleculeName1.length() + MOLECULE_NAME_SEPARATOR.length() + tmpMoleculeName2.length());
                if (tmpMoleculeName1.compareTo(tmpMoleculeName2) < 0) {
                    tmpTupleBuffer.append(tmpMoleculeName1);
                    tmpTupleBuffer.append(MOLECULE_NAME_SEPARATOR);
                    tmpTupleBuffer.append(tmpMoleculeName2);
                } else {
                    tmpTupleBuffer.append(tmpMoleculeName2);
                    tmpTupleBuffer.append(MOLECULE_NAME_SEPARATOR);
                    tmpTupleBuffer.append(tmpMoleculeName1);
                }
                return tmpTupleBuffer.toString();
            default:
                String[] tmpMoleculeNames = aMoleculeNameSet.toArray(new String[0]);
                Arrays.sort(tmpMoleculeNames);
                // Estimated buffer size:
                tmpTupleBuffer = new StringBuilder(2 * aMoleculeNameSet.size() * tmpMoleculeNames[0].length());
                for (int i = 0; i < tmpMoleculeNames.length; i++) {
                    if (i > 0) {
                        tmpTupleBuffer.append(MOLECULE_NAME_SEPARATOR);
                    }
                    tmpTupleBuffer.append(tmpMoleculeNames[i]);
                }
                return tmpTupleBuffer.toString();
        }
    }
    // </editor-fold>
    
}
