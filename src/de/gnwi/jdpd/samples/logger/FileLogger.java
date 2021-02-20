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
package de.gnwi.jdpd.samples.logger;

import de.gnwi.jdpd.utilities.Utils;

/**
 * File logger: Write log information to file when finishing.
 * 
 * @author Achim Zielesny
 */
public class FileLogger extends MemoryLogger {

    // <editor-fold defaultstate="collapsed" desc="Private final class variables">
    /**
     * Path and name of log file
     */
    private final String logFilePathname;
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Constructor">
    /**
     * Constructor
     * 
     * @param aLogFilePathname Log file pathname
     * @param aLogLevels Log levels
     * @throws IllegalArgumentException Thrown if argument is illegal
     */
    public FileLogger(String aLogFilePathname, int[] aLogLevels) throws IllegalArgumentException {
        super(aLogLevels);
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aLogFilePathname == null || aLogFilePathname.isEmpty()) {
            throw new IllegalArgumentException("FileLogger.Constructor: aLogFilePathname is null/empty.");
        }
        // </editor-fold>
        this.logFilePathname = aLogFilePathname;
    }
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public overridden methods">
    /**
     * Finishes logger
     * 
     * @return true: Operation was successful, false: Operation failed
     */
    @Override
    public boolean finish() {
        super.finish();
        if (!Utils.writeStringListToFile(this.logQueue, this.logFilePathname)) {
            return Utils.appendStringListToFile(this.logQueue, this.logFilePathname);
        } else {
            return true;
        }
    }
    // </editor-fold>
    
}
