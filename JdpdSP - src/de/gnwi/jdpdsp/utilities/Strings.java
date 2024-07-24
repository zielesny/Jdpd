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

/**
 * Global strings
 * 
 * @author Achim Zielesny
 */
public final class Strings {

    // <editor-fold defaultstate="collapsed" desc="Time related strings">
    /**
     * Standard timestamp format
     */
    public static final String STANDARD_TIMESTAMP_FORMAT = "yyyy/MM/dd - HH:mm:ss";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Compression related file ending strings">
    /**
     * Ending of GZIP file
     */
    public static String GZIP_FILE_ENDING = ".gz";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Version related strings">
    /**
     * Jdpd version text
     */
    public static final String JDPD_VERSION_TEXT = "Jdpd version " + de.gnwi.jdpd.utilities.Strings.JDPD_VERSION + " (single precision arithmetic)";
    // </editor-fold>
    
}
