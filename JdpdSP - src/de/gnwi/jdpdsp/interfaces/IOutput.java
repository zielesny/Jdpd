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
package de.gnwi.jdpdsp.interfaces;

import de.gnwi.jdpdsp.parameters.RestartInfo;
import de.gnwi.jdpdsp.particlePosition.ParticlePosition;
import de.gnwi.jdpdsp.particlePosition.ParticlePositionPool;
import de.gnwi.jdpdsp.rg.MoleculeRgValue;
import java.util.HashMap;

/**
 * Interface for simulation output
 * 
 * @author Achim Zielesny
 */
public interface IOutput {

    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Sets particle positions at start
     * 
     * @param aParticlePositions Particle positions
     */
    void setStartParticlePositions(ParticlePosition[] aParticlePositions);
    
    /**
     * Sets particle positions of minimization step
     * 
     * @param aMinimizationStep Minimization step
     * @param aParticlePositions Particle positions
     */
    void setMinimizationStepParticlePositions(int aMinimizationStep, ParticlePosition[] aParticlePositions);
    
    /**
     * Sets particle positions after minimization
     * 
     * @param aParticlePositions Particle positions
     */
    void setMinimizedParticlePositions(ParticlePosition[] aParticlePositions);
    
    /**
     * Sets simulation step information
     * 
     * @param aSimulationStep Simulation step
     * @param aTemperature Temperature
     * @param anUpotDpd DPD potential energy
     * @param anUpotBond Bond potential energy
     * @param anUpotElectrostatics Electrostatics potential energy
     * @param anUpotTotal Total potential energy (= aUpotDdpd + aUpotBond + aUpotElectrostatics)
     * @param anUkin Kinetic energy
     * @param anUtotal Total energy
     * @param aSurfaceTensionAlongX Surface tension along x axis
     * @param aSurfaceTensionAlongY Surface tension along y axis
     * @param aSurfaceTensionAlongZ Surface tension along z axis
     * @param aSurfaceTensionNorm Norm (magnitude) of surface tension
     * @param aDpdSurfaceTensionAlongX DPD surface tension along x axis
     * @param aDpdSurfaceTensionAlongY DPD surface tension along y axis
     * @param aDpdSurfaceTensionAlongZ DPD surface tension along z axis
     * @param aDpdSurfaceTensionNorm Norm (magnitude) of DPD surface tension
     * @param aMoleculeRgValues Molecule Rg values (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap Base molecule-particle to nearest-neighbor molecule-particle frequency map (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap Base molecule-particle to nearest-neighbor particle frequency map (may be null)
     * @param aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap Base molecule-particle to nearest-neighbor molecule frequency map (may be null)
     * @param aBaseMoleculeToNearestNeighborMoleculeFrequencyMap Base molecule to nearest-neighbor molecule frequency map (may be null)
     * @param aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap Base molecule to nearest-neighbor molecule-tuple frequency map (may be null)
     * @param aParticlePositions Particle positions
     */
    void setSimulationStepInformation(
        int aSimulationStep, 
        float aTemperature, 
        float anUpotDpd,
        float anUpotBond,
        float anUpotElectrostatics,
        float anUpotTotal, 
        float anUkin, 
        float anUtotal, 
        float aSurfaceTensionAlongX,
        float aSurfaceTensionAlongY,
        float aSurfaceTensionAlongZ,
        float aSurfaceTensionNorm,
        float aDpdSurfaceTensionAlongX,
        float aDpdSurfaceTensionAlongY,
        float aDpdSurfaceTensionAlongZ,
        float aDpdSurfaceTensionNorm,
        MoleculeRgValue[] aMoleculeRgValues,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborMoleculeParticleFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborParticleFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeParticleToNearestNeighborMoleculeFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeToNearestNeighborMoleculeFrequencyMap,
        HashMap<String, HashMap<String, Integer>> aBaseMoleculeToNearestNeighborMoleculeTupleFrequencyMap,
        ParticlePosition[] aParticlePositions
    );
    
    /**
     * Sets restart info
     * 
     * @param aRestartInfo Restart info
     */
    void setRestartInfo(RestartInfo aRestartInfo);
    
    /**
     * Particle position pool
     * 
     * @param aParticlePositionPool Particle position pool
     */
    void setParticlePositionPool(ParticlePositionPool aParticlePositionPool);
    
    /**
     * Finishes output
     * 
     * @return True: Successful finished, false: Otherwise
     */
    boolean finish();
    // </editor-fold>
    
}
