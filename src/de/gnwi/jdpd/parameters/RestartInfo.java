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
package de.gnwi.jdpd.parameters;

import de.gnwi.jdpd.utilities.Constants;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * Restart information
 * 
 * @author Achim Zielesny
 */
public class RestartInfo {
    
    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * Version
     */
    private final static String VERSION_1_0_0 = "Version 1.0.0";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Last simulation time step to which the particle positions and velocities belong
     */
    private final int lastTimeStep;

    /**
     * Number of additional simulation time steps
     */
    private final int additionalTimeStepNumber;
    
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
     * Current x-components of particle velocities
     */
    private final double[] v_x;

    /**
     * Current y-components of particle velocities
     */
    private final double[] v_y;
    
    /**
     * Current z-components of particle velocities
     */
    private final double[] v_z;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * (No checks are performed)
     * 
     * @param aLastTimeStep Last simulation time step to which the particle positions and velocities belong
     * @param aR_x Current x-position of particle in simulation box
     * @param aR_y Current y-position of particle in simulation box
     * @param aR_z Current z-position of particle in simulation box
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     */
    public RestartInfo(
        int aLastTimeStep,
        double[] aR_x,
        double[] aR_y,
        double[] aR_z,
        double[] aV_x,
        double[] aV_y,
        double[] aV_z) {
        this.lastTimeStep = aLastTimeStep;
        this.additionalTimeStepNumber = -1;
        this.r_x = aR_x;
        this.r_y = aR_y;
        this.r_z = aR_z;
        this.v_x = aV_x;
        this.v_y = aV_y;
        this.v_z = aV_z;
    }

    /**
     * Constructor
     * 
     * @param anAdditionalTimeStepNumber Number of additional simulation time steps
     * @param aRestartInfoFilePathname Pathname of file written with this.writeToFile()
     * @throws IllegalArgumentException Thrown if an argument is illegal.
     */
    public RestartInfo(
        int anAdditionalTimeStepNumber,
        String aRestartInfoFilePathname
        ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anAdditionalTimeStepNumber < 1) {
            throw new IllegalArgumentException("RestartInfo.Constructor: anAdditionalTimeStepNumber < 1.");
        }
        if (aRestartInfoFilePathname == null || aRestartInfoFilePathname.isEmpty()) {
            throw new IllegalArgumentException("RestartInfo.Constructor: aRestartInfoFilePathname is null/empty.");
        }
        if (!(new File(aRestartInfoFilePathname)).isFile()) {
            throw new IllegalArgumentException("RestartInfo.Constructor: aRestartInfoFilePathname does not exist.");
        }
        // </editor-fold>
        this.additionalTimeStepNumber = anAdditionalTimeStepNumber;
        BufferedReader tmpBufferedReader = null;
        try {
            String tmpLine;
            tmpBufferedReader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(aRestartInfoFilePathname), Constants.BUFFER_SIZE)));
            tmpLine = tmpBufferedReader.readLine();
            if (tmpLine.equals(VERSION_1_0_0)) {
                // <editor-fold defaultstate="collapsed" desc="Version 1.0.0">
                this.lastTimeStep = Integer.valueOf(tmpBufferedReader.readLine());
                int tmpParticleNumber = Integer.valueOf(tmpBufferedReader.readLine());
                this.r_x = new double[tmpParticleNumber];
                this.r_y = new double[tmpParticleNumber];
                this.r_z = new double[tmpParticleNumber];
                this.v_x = new double[tmpParticleNumber];
                this.v_y = new double[tmpParticleNumber];
                this.v_z = new double[tmpParticleNumber];
                for (int i = 0; i < tmpParticleNumber; i++) {
                    this.r_x[i] = Double.valueOf(tmpBufferedReader.readLine());
                    this.r_y[i] = Double.valueOf(tmpBufferedReader.readLine());
                    this.r_z[i] = Double.valueOf(tmpBufferedReader.readLine());
                    this.v_x[i] = Double.valueOf(tmpBufferedReader.readLine());
                    this.v_y[i] = Double.valueOf(tmpBufferedReader.readLine());
                    this.v_z[i] = Double.valueOf(tmpBufferedReader.readLine());
                }
                // </editor-fold>
            } else {
                throw new IllegalArgumentException("RestartInfo.Constructor: Unknown version.");
            }
        } catch (Exception anException) {
            throw new IllegalArgumentException("RestartInfo.Constructor: Can not read file.");
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException e) {
                    throw new IllegalArgumentException("RestartInfo.Constructor: Can not close file.");
                }
            }
        }
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    /**
     * Write restart info to file
     * 
     * @param aRestartInfoFilePathname Pathname of restart info file
     * @return True: Operation successful, false: Otherwise
     */
    public boolean writeToFile(String aRestartInfoFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aRestartInfoFilePathname == null || aRestartInfoFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        try (PrintWriter tmpPrintWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(aRestartInfoFilePathname), Constants.BUFFER_SIZE))) {
            this.writeOutput(tmpPrintWriter);
        } catch (Exception e) {
            return false;
        }
        return true;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public properties (get)">
    /**
     * Last simulation time step to which the particle positions and velocities belong
     * 
     * @return Last simulation time step to which the particle positions and velocities belong
     */
    public int getLastTimeStep () {
        return this.lastTimeStep;
    }

    /**
     * Number of additional simulation time steps
     * 
     * @return Number of additional simulation time steps
     */
    public int getAdditionalTimeStepNumber() {
        return this.additionalTimeStepNumber;
    }
    
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
     * Current x-components of particle velocities
     * 
     * @return Current x-components of particle velocities
     */
    public double[] getV_x() {
        return this.v_x;
    };

    /**
     * Current y-components of particle velocities
     * 
     * @return Current y-components of particle velocities
     */
    public double[] getV_y() {
        return this.v_y;
    };

    /**
     * Current z-components of particle velocities
     * 
     * @return Current z-components of particle velocities
     */
    public double[] getV_z() {
        return this.v_z;
    };
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Writes output (see code).
     * NOTE: No checks are performed!
     * 
     * @param aPrintWriter Print writer to write to
     */
    private void writeOutput(PrintWriter aPrintWriter) {
        aPrintWriter.println(VERSION_1_0_0);
        aPrintWriter.println(String.valueOf(this.lastTimeStep));
        aPrintWriter.println(String.valueOf(this.r_x.length));
        for (int i = 0; i < this.r_x.length; i++) {
            aPrintWriter.println(String.valueOf(this.r_x[i]));
            aPrintWriter.println(String.valueOf(this.r_y[i]));
            aPrintWriter.println(String.valueOf(this.r_z[i]));
            aPrintWriter.println(String.valueOf(this.v_x[i]));
            aPrintWriter.println(String.valueOf(this.v_y[i]));
            aPrintWriter.println(String.valueOf(this.v_z[i]));
        }
    }
    // </editor-fold>
    
}
