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

/**
 * Global file output strings
 * 
 * @author Achim Zielesny
 */
public final class FileOutputStrings {

    /**
     * Particle positions start filename prefix
     */
    public final static String PARTICLE_POSITIONS_START_FILENAME_PREFIX = "PPstart";

    /**
     * Format for particle positions minimization step file prefix
     */
    public final static String PARTICLE_POSITIONS_MINIMIZATION_STEP_FILE_PREFIX_FORMAT = "%sPPmin%s";

    /**
     * Minimised particle positions filename prefix 
     */
    public final static String PARTICLE_POSITIONS_MINIMIZED_FILENAME_PREFIX = "PPmin";

    /**
     * Final particle positions filename prefix (after simulation)
     */
    public final static String PARTICLE_POSITIONS_FINAL_FILENAME_PREFIX = "PPfinal";

    /**
     * Particle positions simulation step file prefix
     */
    public final static String PARTICLE_POSITIONS_SIMULATION_STEP_FILE_PREFIX = "PP";
    
    /**
     * Format for particle positions simulation step file prefix
     */
    public final static String PARTICLE_POSITIONS_SIMULATION_STEP_FILE_PREFIX_FORMAT = "%sPP%s";

    /**
     * Ending of text file
     */
    public final static String TEXT_FILE_ENDING = ".txt";
    
    /**
     * Property file names
     */
    public final static String SIMULATION_STEP_FILENAME = "SimStep" + TEXT_FILE_ENDING;
    public final static String TEMPERATURE_FILENAME = "T" + TEXT_FILE_ENDING;
    public final static String U_POT_DPD_FILENAME = "UpotDpd" + TEXT_FILE_ENDING;
    public final static String U_POT_BOND_FILENAME = "UpotBond" + TEXT_FILE_ENDING;
    public final static String U_POT_ELECTROSTATICS_FILENAME = "UpotElectrostatics" + TEXT_FILE_ENDING;
    public final static String U_POT_TOTAL_FILENAME = "UpotTotal" + TEXT_FILE_ENDING;
    public final static String U_KIN_FILENAME = "Ukin" + TEXT_FILE_ENDING;
    public final static String U_TOTAL_FILENAME = "Utotal" + TEXT_FILE_ENDING;
    public final static String SURFACE_TENSION_X_FILENAME = "SurfaceTensionX" + TEXT_FILE_ENDING;
    public final static String SURFACE_TENSION_Y_FILENAME = "SurfaceTensionY" + TEXT_FILE_ENDING;
    public final static String SURFACE_TENSION_Z_FILENAME = "SurfaceTensionZ" + TEXT_FILE_ENDING;
    public final static String SURFACE_TENSION_NORM_FILENAME = "SurfaceTensionNorm" + TEXT_FILE_ENDING;
    public final static String DPD_SURFACE_TENSION_X_FILENAME = "DpdSurfaceTensionX" + TEXT_FILE_ENDING;
    public final static String DPD_SURFACE_TENSION_Y_FILENAME = "DpdSurfaceTensionY" + TEXT_FILE_ENDING;
    public final static String DPD_SURFACE_TENSION_Z_FILENAME = "DpdSurfaceTensionZ" + TEXT_FILE_ENDING;
    public final static String DPD_SURFACE_TENSION_NORM_FILENAME = "DpdSurfaceTensionNorm" + TEXT_FILE_ENDING;

    /**
     * Prefix of radius-of-gyration file
     */
    public final static String RADIUS_OF_GYRATION_FILENAME_PREFIX = "Rg_";

    /**
     * Nearest-neighbor file names
     */
    public final static String MP_TO_NN_MP_FILENAME = "MP_MP" + TEXT_FILE_ENDING;
    public final static String MP_TO_NN_P_FILENAME = "MP_P" + TEXT_FILE_ENDING;
    public final static String MP_TO_NN_M_FILENAME = "MP_M" + TEXT_FILE_ENDING;
    public final static String M_TO_NN_M_FILENAME = "M_M" + TEXT_FILE_ENDING;
    public final static String M_TO_NN_M_TUPLE_FILENAME = "M_M_TUPLE" + TEXT_FILE_ENDING;
            
    /**
     * Restart info file name prefix
     */
    public final static String RESTART_INFO_FILENAME_PREFIX = "RestartInfo";
    
}
