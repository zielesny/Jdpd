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
package de.gnwi.jdpd.utilities;

/**
 * Global constants
 * 
 * @author Achim Zielesny
 */
public final class Constants {

    // <editor-fold defaultstate="collapsed" desc="Buffer related definitions">
    /**
     * Buffer size (64 kByte = 65536, 256 kByte = 262144, 512 kByte = 524288, 1
     * MByte = 1048576 Byte)
     */
    public static final int BUFFER_SIZE = 65536;
    // </editor-fold>
    
    // <editor-fold defaultstate="collapsed" desc="Numeric definitions">
    /**
     * Maximum number of time steps for a simulation
     */
    public static final int MAXIMUM_NUMBER_OF_TIME_STEPS = 1000000000;

    /**
     * Tiny threshold (a little smaller than the significant number of digits 
     * of single (!) precision)
     */
    public static final double TINY_THRESHOLD = 1E-6;
    // </editor-fold>
    
}
