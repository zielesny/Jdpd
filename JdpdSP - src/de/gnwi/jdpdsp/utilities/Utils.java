/**
 * JdpdSP - Molecular Fragment Dissipative Particle Dynamics (DPD) Simulation
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
package de.gnwi.jdpdsp.utilities;

import de.gnwi.jdpdsp.accumulators.ForceAccumulator;
import de.gnwi.jdpdsp.accumulators.PotentialAccumulator;
import de.gnwi.jdpdsp.movement.MoleculeFixationInfo;
import de.gnwi.jdpdsp.movement.MoleculeVelocityFixationInfo;
import de.gnwi.jdpdsp.interfaces.ILogger;
import de.gnwi.jdpdsp.interfaces.IOutput;
import de.gnwi.jdpdsp.movement.MoleculeAccelerationInfo;
import de.gnwi.jdpdsp.movement.MoleculeBoundaryInfo;
import de.gnwi.jdpdsp.movement.MoleculeSphereInfo;
import de.gnwi.jdpdsp.parameters.ChemicalSystemDescription;
import de.gnwi.jdpdsp.particlePosition.ParticlePositionPool;
import de.gnwi.jdpdsp.particlePosition.ParticlePosition;
import de.gnwi.jdpdsp.parameters.ParticleArrays;
import de.gnwi.jdpdsp.parameters.MoleculeTypes;
import de.gnwi.jdpdsp.parameters.Parameters;
import de.gnwi.jdpdsp.parameters.SimulationDescription;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.math3.util.FastMath;

/**
 * Thread-safe static utility methods
 * 
 * @author Achim Zielesny
 */
public final class Utils {

    // <editor-fold defaultstate="collapsed" desc="Private static final class variables">
    /**
     * General separator string
     */
    private static final String GENERAL_SEPARATOR = "|";

    /**
     * Time difference format strings
     */
    private static final String FORMAT_TIME_DIFFERENCE_SECONDS = "%s s";
    private static final String FORMAT_TIME_DIFFERENCE_MINUTES = "%s min %s s";
    private static final String FORMAT_TIME_DIFFERENCE_HOURS = "%s h %s min %s s";
    private static final String FORMAT_TIME_DIFFERENCE_DAYS = "%s d %s h %s min %s s";

    /**
     * Time factors
     */
    private static final long DAY_FACTOR = 1000L * 60L * 60L * 24L;
    private static final long HOUR_FACTOR = 1000L * 60L * 60L;
    private static final long MINUTES_FACTOR = 1000L * 60L;
    private static final long SECONDS_FACTOR = 1000L;

    /**
     * Numeric constants
     */
    private static final float THREE = 3.0f;
    private static final float A_HALF = 0.5f;
    private static final float ONE_A_HALF = 1.5f;
    private final static float ONE = 1.0f;
    private final static float TWO = 2.0f;
    private final static float SEVEN = 7.0f;
    private final static float FOURTY_TWO = 42.0f;
    private final static float C7_4 = 7.0f / 4.0f;
    private final static float C21_4 = 21.0f / 4.0f;
    private final static float C5_2 = 5.0f / 2.0f;
    private final static float C35_2 = 35.0f / 2.0f;
    private final static float C105_4 = 105.0f / 4.0f;

    /**
     * Separator string for molecule name and particle token
     */
    private static final String MOLECULE_NAME_PARTICLE_TOKEN_SEPARATOR = "-";
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Public methods">
    // <editor-fold defaultstate="collapsed" desc="File related methods">
    /**
     * Reads string list from file
     *
     * @param aSourceFilePathname Full pathname of source (may be null/empty
     * then null is returned)
     * @param aCommentLinePrefix Prefix of comment line (may be null/empty)
     * @return String list or null if string array could not be read
     */
    public static LinkedList<String> readStringListFromFile(String aSourceFilePathname, String aCommentLinePrefix) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty() || !(new File(aSourceFilePathname)).isFile()) {
            return null;
        }
        // </editor-fold>
        BufferedReader tmpBufferedReader = null;
        try {
            if (!(new File(aSourceFilePathname)).isFile()) {
                return null;
            }
            FileReader tmpFileReader = new FileReader(aSourceFilePathname);
            tmpBufferedReader = new BufferedReader(tmpFileReader, Constants.BUFFER_SIZE);
            LinkedList<String> tmpLinkedList = new LinkedList<>();
            String tmpLine;
            if (aCommentLinePrefix == null || aCommentLinePrefix.isEmpty()) {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    tmpLinkedList.add(tmpLine);
                }
            } else {
                while ((tmpLine = tmpBufferedReader.readLine()) != null) {
                    if (!tmpLine.startsWith(aCommentLinePrefix)) {
                        tmpLinkedList.add(tmpLine);
                    }
                }
            }
            return tmpLinkedList;
        } catch (Exception e) {
            return null;
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException anException) {
                    return null;
                }
            }
        }
    }
    
    /**
     * Writes string list to file.
     * NOTE: If destination file pathname already exists nothing is done.
     *
     * @param aStringList String list (may be null or empty then false is
     * returned)
     * @param aDestinationFilePathname Full pathname of destination (may be null
     * then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean writeStringListToFile(Iterable<String> aStringList, String aDestinationFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aStringList == null || aDestinationFilePathname == null || aDestinationFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        PrintWriter tmpPrintWriter = null;
        try {
            if ((new File(aDestinationFilePathname)).isFile()) {
                return false;
            }
            FileWriter tmpFileWriter = new FileWriter(aDestinationFilePathname);
            BufferedWriter tmpBufferedWriter = new BufferedWriter(tmpFileWriter, Constants.BUFFER_SIZE);
            tmpPrintWriter = new PrintWriter(tmpBufferedWriter);
            for (String tmpSingleString : aStringList) {
                if (tmpSingleString != null) {
                    tmpPrintWriter.println(tmpSingleString);
                } else {
                    // If string is null then write empty string to guarantee
                    // same size of list in read operation
                    tmpPrintWriter.println("");
                }
            }
            tmpPrintWriter.flush();
            return true;
        } catch (Exception anException) {
            return false;
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.close();
            }
        }
    }
    
    /**
     * Append string list to existing file.
     * NOTE: If destination file pathname does not exist nothing is done.
     *
     * @param aStringList String list (may be null or empty then false is
     * returned)
     * @param aDestinationFilePathname Full pathname of destination (may be null
     * then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean appendStringListToFile(Iterable<String> aStringList, String aDestinationFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aStringList == null || aDestinationFilePathname == null || aDestinationFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        PrintWriter tmpPrintWriter = null;
        try {
            if (!(new File(aDestinationFilePathname)).isFile()) {
                return false;
            }
            // Parameter true: Append
            FileWriter tmpFileWriter = new FileWriter(aDestinationFilePathname, true);
            BufferedWriter tmpBufferedWriter = new BufferedWriter(tmpFileWriter, Constants.BUFFER_SIZE);
            tmpPrintWriter = new PrintWriter(tmpBufferedWriter);
            for (String tmpSingleString : aStringList) {
                if (tmpSingleString != null) {
                    tmpPrintWriter.println(tmpSingleString);
                } else {
                    // If string is null then write empty string to guarantee
                    // same size of list in read operation
                    tmpPrintWriter.println("");
                }
            }
            tmpPrintWriter.flush();
            return true;
        } catch (Exception anException) {
            return false;
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.close();
            }
        }
    }

    /**
     * Reads string array from file. NOTE: First line must be size of string
     * array.
     *
     * @param aSourceFilePathname Full pathname of source (may be null then null
     * is returned)
     * @return String array or null if string array could not be read
     */
    public static String[] readStringArrayFromFile(String aSourceFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty()) {
            return null;
        }
        // </editor-fold>
        BufferedReader tmpBufferedReader = null;
        try {
            if (!(new File(aSourceFilePathname)).isFile()) {
                return null;
            }
            FileReader tmpFileReader = new FileReader(aSourceFilePathname);
            tmpBufferedReader = new BufferedReader(tmpFileReader, Constants.BUFFER_SIZE);
            int tmpSize = Integer.valueOf(tmpBufferedReader.readLine());
            String[] tmpStringArray = new String[tmpSize];
            for (int i = 0; i < tmpSize; i++) {
                tmpStringArray[i] = tmpBufferedReader.readLine();
            }
            return tmpStringArray;
        } catch (Exception anException) {
            return null;
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException e) {
                    return null;
                }
            }
        }
    }

    /**
     * Writes string array to file. First line is the number of lines. NOTE: If
     * destination pathname already exists nothing is done.
     *
     * @param aStringArray String array (may be null then false is returned)
     * @param aDestinationFilePathname Full pathname of destination (may be null
     * then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean writeStringArrayToFile(String[] aStringArray, String aDestinationFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aStringArray == null || aStringArray.length == 0 || aDestinationFilePathname == null || aDestinationFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        PrintWriter tmpPrintWriter = null;
        try {
            if ((new File(aDestinationFilePathname)).isFile()) {
                return false;
            }
            FileWriter tmpFileWriter = new FileWriter(aDestinationFilePathname);
            BufferedWriter tmpBufferedWriter = new BufferedWriter(tmpFileWriter, Constants.BUFFER_SIZE);
            tmpPrintWriter = new PrintWriter(tmpBufferedWriter);
            // Write number of lines
            tmpPrintWriter.println(String.valueOf(aStringArray.length));
            for (int i = 0; i < aStringArray.length; i++) {
                if (aStringArray[i] != null) {
                    tmpPrintWriter.println(aStringArray[i]);
                } else {
                    // If string is null then write empty string to guarantee
                    // same dimension of string array in read operation
                    tmpPrintWriter.println("");
                }
            }
            tmpPrintWriter.flush();
            return true;
        } catch (Exception anException) {
            return false;
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.close();
            }
        }
    }

    /**
     * Writes aBaseToNearestNeighborStepFrequencyMap to file.
     * NOTE: Bases and nearest neighbors are sorted ascending in file.
     * NOTE: If destination pathname already exists nothing is done.
     *
     * @param aBaseToNearestNeighborStepFrequencyMap Base to nearest-neighbor 
     * step-frequency map (may be null then false is returned)
     * @param aDestinationFilePathname Full pathname of destination (may be null
     * then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean writeBaseToNearestNeighborStepFrequencyMapToFile(
        HashMap<String, HashMap<String, LinkedList<String>>> aBaseToNearestNeighborStepFrequencyMap, 
        String aDestinationFilePathname
    ) 
    {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aBaseToNearestNeighborStepFrequencyMap == null 
            || aBaseToNearestNeighborStepFrequencyMap.isEmpty() 
            || aDestinationFilePathname == null 
            || aDestinationFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        PrintWriter tmpPrintWriter = null;
        try {
            if ((new File(aDestinationFilePathname)).isFile()) {
                return false;
            }
            FileWriter tmpFileWriter = new FileWriter(aDestinationFilePathname);
            BufferedWriter tmpBufferedWriter = new BufferedWriter(tmpFileWriter, Constants.BUFFER_SIZE);
            tmpPrintWriter = new PrintWriter(tmpBufferedWriter);
            // Write version
            tmpPrintWriter.println(de.gnwi.jdpd.utilities.Strings.INPUT_OUTPUT_VERSION_1_0_0);
            // Sort bases
            String[] tmpBaseArray = aBaseToNearestNeighborStepFrequencyMap.keySet().toArray(new String[0]);
            Arrays.sort(tmpBaseArray);
            // Write number of bases
            tmpPrintWriter.println(String.valueOf(tmpBaseArray.length));
            for (String tmpBase : tmpBaseArray) {
                // Write base
                tmpPrintWriter.println(tmpBase);
                HashMap<String, LinkedList<String>> tmpNearestNeighborToStepFrequencyMap = 
                    aBaseToNearestNeighborStepFrequencyMap.get(tmpBase);
                // Sort nearest neighbors
                String[] tmpNearestNeighborArray = tmpNearestNeighborToStepFrequencyMap.keySet().toArray(new String[0]);
                Arrays.sort(tmpNearestNeighborArray);
                // Write number of nearest neighbors
                tmpPrintWriter.println(String.valueOf(tmpNearestNeighborArray.length));
                for (String tmpNearestNeighbor : tmpNearestNeighborArray) {
                    // Write nearest neighbor
                    tmpPrintWriter.println(tmpNearestNeighbor);
                    LinkedList<String> tmpStepFrequencyList = tmpNearestNeighborToStepFrequencyMap.get(tmpNearestNeighbor);
                    // Write number of list elements
                    // NOTE: First list element is version
                    tmpPrintWriter.println(String.valueOf(tmpStepFrequencyList.size()));
                    for (String tmpListElement : tmpStepFrequencyList) {
                        // Write list element
                        tmpPrintWriter.println(tmpListElement);
                    }
                }
            }
            tmpPrintWriter.flush();
            return true;
        } catch (Exception anException) {
            return false;
        } finally {
            if (tmpPrintWriter != null) {
                tmpPrintWriter.close();
            }
        }
    }

    /**
     * Reads baseToNearestNeighborStepFrequencyMap from file which was written
     * with method writeBaseToNearestNeighborStepFrequencyMapToFile().
     *
     * @param aSourceFilePathname Full pathname of source (may be null then null
     * is returned)
     * @return Hash map baseToNearestNeighborStepFrequencyMap or null if 
     * information could not be read
     */
    public static HashMap<String, HashMap<String, LinkedList<String>>> readBaseToNearestNeighborStepFrequencyMapFromFile(String aSourceFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty()) {
            return null;
        }
        // </editor-fold>
        BufferedReader tmpBufferedReader = null;
        try {
            if (!(new File(aSourceFilePathname)).isFile()) {
                return null;
            }
            FileReader tmpFileReader = new FileReader(aSourceFilePathname);
            tmpBufferedReader = new BufferedReader(tmpFileReader, Constants.BUFFER_SIZE);
            // Read version
            String tmpVersion = tmpBufferedReader.readLine();
            if (!tmpVersion.equals(de.gnwi.jdpd.utilities.Strings.INPUT_OUTPUT_VERSION_1_0_0)) {
                return null;
            }
            // Read number of bases
            int tmpNumberOfBases = Integer.valueOf(tmpBufferedReader.readLine());
            HashMap<String, HashMap<String, LinkedList<String>>> tmpBaseToNearestNeighborStepFrequencyMap = 
                new HashMap<>(tmpNumberOfBases);
            for (int tmpBaseCounter = 0; tmpBaseCounter < tmpNumberOfBases; tmpBaseCounter++) {
                // Read base
                String tmpBase = tmpBufferedReader.readLine();
                // Read number of nearest neighbors
                int tmpNumberOfNearestNeighbors = Integer.valueOf(tmpBufferedReader.readLine());
                HashMap<String, LinkedList<String>> tmpNearestNeighborToStepFrequencyMap = new HashMap<>(tmpNumberOfNearestNeighbors);
                for (int tmpNearestNeighborCounter = 0; tmpNearestNeighborCounter < tmpNumberOfNearestNeighbors; tmpNearestNeighborCounter++) {
                    // Read nearest neighbor
                    String tmpNearestNeighbor = tmpBufferedReader.readLine();
                    // Read number of list elements
                    int tmpNumberOfListElements = Integer.valueOf(tmpBufferedReader.readLine());
                    LinkedList<String> tmpStepFrequencyList = new LinkedList<>();
                    for (int tmpListElementCounter = 0; tmpListElementCounter < tmpNumberOfListElements; tmpListElementCounter++) {
                        // Read list element
                        tmpStepFrequencyList.add(tmpBufferedReader.readLine());
                    }
                    tmpNearestNeighborToStepFrequencyMap.put(tmpNearestNeighbor, tmpStepFrequencyList);
                }
                tmpBaseToNearestNeighborStepFrequencyMap.put(tmpBase, tmpNearestNeighborToStepFrequencyMap);
            }
            return tmpBaseToNearestNeighborStepFrequencyMap;
        } catch (Exception anException) {
            return null;
        } finally {
            if (tmpBufferedReader != null) {
                try {
                    tmpBufferedReader.close();
                } catch (IOException e) {
                    return null;
                }
            }
        }
    }

    /**
     * Deletes defined files and directories
     *
     * @param aFilesAndDirectoriesForDeletion Array of defined files and directories
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean deleteDefinedFilesAndDirectories(File[] aFilesAndDirectoriesForDeletion) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aFilesAndDirectoriesForDeletion == null || aFilesAndDirectoriesForDeletion.length == 0) {
            return true;
        }
        // </editor-fold>
        try {
            for (File file : aFilesAndDirectoriesForDeletion) {
                if (file.isDirectory()) {
                    if (!Utils.deleteDirectory(file)) {
                        return false;
                    }
                } else {
                    if (!file.delete()) {
                        return false;
                    }
                }
            }
            return true;
        } catch (Exception anException) {
            return false;
        }
    }

    /**
     * Deletes directory and all sub directories
     *
     * @param aDirectory Directory object to be deleted (may be null then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean deleteDirectory(File aDirectory) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aDirectory == null) {
            return false;
        }
        if (!aDirectory.isDirectory()) {
            // Directory does not exist, i.e. it is deleted: Return true
            return true;
        }
        // </editor-fold>
        try {
            if (!Utils.deleteDefinedFilesAndDirectories(aDirectory.listFiles())) {
                return false;
            }
            return aDirectory.delete();
        } catch (Exception anException) {
            return false;
        }
    }

    /**
     * Deletes directory and all sub directories
     *
     * @param aDirectoryPath Path of directory to be deleted (may be null or
     * empty then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean deleteDirectory(String aDirectoryPath) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aDirectoryPath == null || aDirectoryPath.isEmpty()) {
            return false;
        }
        if (!(new File(aDirectoryPath)).isDirectory()) {
            // Directory does not exist, i.e. it is deleted: Return true
            return true;
        }
        // </editor-fold>
        try {
            return Utils.deleteDirectory(new File(aDirectoryPath));
        } catch (Exception anException) {
            return false;
        }
    }

    /**
     * Clears directory from all files and all sub directories
     *
     * @param aDirectory Directory object to be cleared (may be null then false
     * is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean clearDirectory(File aDirectory) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aDirectory == null) {
            return false;
        }
        if (!aDirectory.isDirectory()) {
            // Directory does not exist
            return false;
        }
        // </editor-fold>
        return Utils.deleteDefinedFilesAndDirectories(aDirectory.listFiles());
    }

    /**
     * Clears directory from all files and all sub directories
     *
     * @param aDirectoryPath Directory path to be cleared (may be null then false
     * is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean clearDirectory(String aDirectoryPath) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aDirectoryPath == null || aDirectoryPath.isEmpty()) {
            return false;
        }
        if (!(new File(aDirectoryPath)).isDirectory()) {
            // Directory does not exist, i.e. it is deleted: Return true
            return true;
        }
        // </editor-fold>
        return Utils.clearDirectory(new File(aDirectoryPath));
    }
    
    /**
     * Deletes single file
     *
     * @param aFilePathname Full pathname of file to be deleted (may be null
     * then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean deleteSingleFile(String aFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aFilePathname == null || aFilePathname.isEmpty()) {
            return false;
        }
        // </editor-fold>
        try {
            File tmpFile = new File(aFilePathname);
            if (!tmpFile.isFile()) {
                return true;
            } else {
                return tmpFile.delete();
            }
        } catch (Exception e) {
            return false;
        }
    }
    
    /**
     * Copies specified file to specified destination. NOTE: File will not be
     * copied if specified destination pathname already exists. Then true is
     * returned.
     *
     * @param aSourceFilePathname Full pathname of source (may be null or empty
     * then false is returned)
     * @param aDestinationFilePathname Full pathname of destination (may be null
     * or empty then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean copySingleFile(String aSourceFilePathname, String aDestinationFilePathname) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty() || aDestinationFilePathname == null || aDestinationFilePathname.isEmpty()) {
            return false;
        }
        if (aSourceFilePathname.equals(aDestinationFilePathname)) {
            // File is already in desired destination: Return true
            return true;
        }
        if (!(new File(aSourceFilePathname)).isFile()) {
            return false;
        }
        if ((new File(aDestinationFilePathname)).isFile()) {
            return false;
        }
        // </editor-fold>
        FileInputStream tmpFileInputStream = null;
        FileOutputStream tmpFileOutputStream = null;
        try {
            tmpFileInputStream = new FileInputStream(aSourceFilePathname);
            tmpFileOutputStream = new FileOutputStream(aDestinationFilePathname);
            int tmpLength;
            byte[] tmpBuffer = new byte[Constants.BUFFER_SIZE];
            while ((tmpLength = tmpFileInputStream.read(tmpBuffer)) != -1) {
                tmpFileOutputStream.write(tmpBuffer, 0, tmpLength);
            }
            return true;
        } catch (Exception e) {
            return false;
        } finally {
            try {
                if (tmpFileInputStream != null) {
                    tmpFileInputStream.close();
                }
                if (tmpFileOutputStream != null) {
                    tmpFileOutputStream.close();
                }
            } catch (Exception e) {
                return false;
            }
        }
    }
        
    /**
     * Copies file defined in aSourceFilePathname to temporary directory of OS 
     * and renames the file to aNewFilename.
     * If renamed file already exists in temporary directory of OS it is 
     * overwritten.
     * 
     * @param aSourceFilePathname Full pathname of source (may be null or empty
     * then false is returned)
     * @param aNewFilename New file name (may be null or empty then false is 
     * returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean copySingleFileToTemporaryDirectory(String aSourceFilePathname, String aNewFilename) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSourceFilePathname == null || aSourceFilePathname.isEmpty()) {
            return false;
        }
        if (aNewFilename == null || aNewFilename.isEmpty()) {
            return false;
        }
        // </editor-fold>
        try {
            File tmpSourceFile = new File(aSourceFilePathname);
            if (!tmpSourceFile.isFile()) {
                return false;
            }
            String tmpTemporaryDirectory = System.getProperty("java.io.tmpdir");
            String tmpDestinationFilePathname = null;
            if (tmpTemporaryDirectory.endsWith(String.valueOf(File.separatorChar))) {
                tmpDestinationFilePathname = tmpTemporaryDirectory + aNewFilename;
            } else {
                tmpDestinationFilePathname = tmpTemporaryDirectory + File.separatorChar + aNewFilename;
            }
            if ((new File(tmpDestinationFilePathname)).exists()) {
                (new File(tmpDestinationFilePathname)).delete();
            }
            return Utils.copySingleFile(aSourceFilePathname, tmpDestinationFilePathname);
        } catch (Exception e) {
            return false;
        }
    }
    
    /**
     * Moves file from temporary directory of OS to destination.
     * 
     * @param aFilename Name of file to be moved (may be null or empty then 
     * false is returned)
     * @param aDestinationDirectoryPath Full path of destination directory the
     * file is moved to (may be null then false is returned)
     * @return true: Operation was successful, false: Operation failed
     */
    public static boolean moveSingleFileFromTemporaryDirectory(String aFilename, String aDestinationDirectoryPath) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aFilename == null || aFilename.isEmpty()) {
            return false;
        }
        if (aDestinationDirectoryPath == null || aDestinationDirectoryPath.isEmpty()) {
            return false;
        }
        if (!(new File(aDestinationDirectoryPath)).isDirectory()) {
            return false;
        }
        // </editor-fold>
        try {
            String tmpTemporaryDirectory = System.getProperty("java.io.tmpdir");
            String tmpSourceFilePathname = tmpTemporaryDirectory + File.separatorChar + aFilename;
            if (!(new File(tmpSourceFilePathname)).isFile()) {
                return false;
            }
            String tmpDestinationFilePathname = aDestinationDirectoryPath + File.separatorChar + aFilename;
            if (!Utils.copySingleFile(tmpSourceFilePathname, tmpDestinationFilePathname)) {
                return false;
            }
            if (!(new File(tmpSourceFilePathname)).delete()) {
                return false;
            }
            return true;
        } catch (Exception e) {
            return false;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Exception related methods">
    /**
     * Appends stack trace of exception to string list
     * 
     * @param anException Exception
     * @param aList List
     */
    public static void addStacktrace(Exception anException, List<String> aList) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anException == null) {
            return;
        }
        if (aList == null) {
            return;
        }
        // </editor-fold>
        StringWriter tmpStringWriter = new StringWriter();
        anException.printStackTrace(new PrintWriter(tmpStringWriter));
        String tmpStackTrace = tmpStringWriter.toString();
        String[] tmpStackTraceLines = RegexPatterns.NEW_LINE_PATTERN.split(tmpStackTrace);
        aList.addAll(Arrays.asList(tmpStackTraceLines));
    }
    
    /**
     * Returns stack trace of exception
     * 
     * @param anException Exception
     * @return Stack trace or null if none could be created
     */
    public static String getStacktrace(Exception anException) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (anException == null) {
            return null;
        }
        // </editor-fold>
        StringWriter tmpStringWriter = new StringWriter();
        anException.printStackTrace(new PrintWriter(tmpStringWriter));
        return tmpStringWriter.toString();
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="String related methods">
    /**
     * Splits a string into a string array after one or more whitespace
     * characters.
     *
     * @param aString String to split
     * @return A string array containing the trimmed split strings or null if
     * aString was null or empty
     */
    public static String[] splitAndTrim(String aString) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aString == null || aString.isEmpty()) {
            return null;
        }
        // </editor-fold>
        String[] tmpStringItems = RegexPatterns.WHITESPACE_PATTERN.split(aString);
        if (tmpStringItems == null || tmpStringItems.length == 0) {
            return null;
        }
        // NOTE: Trim-operation is not necessary since all whitespace characters are removed
        return tmpStringItems;
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Time related methods">
    /**
     * Returns current timestamp in format Strings.STANDARD_TIMESTAMP_FORMAT
     *
     * @return Current Timestamp according in format 
     * Strings.STANDARD_TIMESTAMP_FORMAT or null if timestamp can not be 
     * created.
     */
    public static String getTimestamp() {
        SimpleDateFormat tmpSimpleDateFormat;
        try {
            tmpSimpleDateFormat = new SimpleDateFormat(Strings.STANDARD_TIMESTAMP_FORMAT);
        } catch (NullPointerException | IllegalArgumentException e) {
            return null;
        }
        return tmpSimpleDateFormat.format((new GregorianCalendar()).getTime());
    }

    /**
     * Returns formatted time period string
     *
     * @param aTimePeriodInMillis Time period in milliseconds
     * @return Time period string in defined format or null if none could be
     * generated
     */
    public static String getFormattedTimePeriodString(long aTimePeriodInMillis) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aTimePeriodInMillis < 0) {
            return null;
        }
        // </editor-fold>
        try {
            int tmpDays = (int) FastMath.floor(aTimePeriodInMillis / DAY_FACTOR);
            long tmpHourRemainder = aTimePeriodInMillis - tmpDays * DAY_FACTOR;
            int tmpHours = (int) FastMath.floor(tmpHourRemainder / HOUR_FACTOR);
            long tmpMinRemainder = tmpHourRemainder - tmpHours * HOUR_FACTOR;
            int tmpMinutes = (int) FastMath.floor(tmpMinRemainder / MINUTES_FACTOR);
            long tmpSecondsRemainder = tmpMinRemainder - tmpMinutes * MINUTES_FACTOR;
            int tmpSeconds = (int) FastMath.floor(tmpSecondsRemainder / SECONDS_FACTOR);
            if (tmpDays == 0) {
                if (tmpHours == 0) {
                    if (tmpMinutes == 0) {
                        return String.format(FORMAT_TIME_DIFFERENCE_SECONDS, Integer.toString(tmpSeconds));
                    } else {
                        return String.format(FORMAT_TIME_DIFFERENCE_MINUTES, Integer.toString(tmpMinutes), Integer.toString(tmpSeconds));
                    }
                } else {
                    return String.format(FORMAT_TIME_DIFFERENCE_HOURS, Integer.toString(tmpHours), Integer.toString(tmpMinutes), Integer.toString(tmpSeconds));
                }
            } else {
                return String.format(FORMAT_TIME_DIFFERENCE_DAYS, Integer.toString(tmpDays), Integer.toString(tmpHours), Integer.toString(tmpMinutes), Integer.toString(tmpSeconds));
            }
        } catch (Exception e) {
            return null;
        }
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Array copy methods">
    /**
     * Copies arrays to old arrays
     * (No checks are performed)
     * 
     * @param aX x array
     * @param aY y array
     * @param aZ z array
     * @param aOldX Old x array (may be changed)
     * @param aOldY Old y array (may be changed)
     * @param aOldZ Old z array (may be changed)
     */
    public static void copyToOld(
        float[] aX,
        float[] aY,
        float[] aZ,
        float[] aOldX,
        float[] aOldY,
        float[] aOldZ
    ) {
        System.arraycopy(
            aX, 
            0, 
            aOldX, 
            0, 
            aX.length
        );
        System.arraycopy(
            aY, 
            0, 
            aOldY, 
            0, 
            aY.length
        );
        System.arraycopy(
            aZ, 
            0, 
            aOldZ, 
            0, 
            aZ.length
        );
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Calculation methods">
    /**
     * Adds b to a
     * (No checks are performed)
     * 
     * @param a_x x-component of a (may be changed)
     * @param a_y y-component of a (may be changed)
     * @param a_z z-component of a (may be changed)
     * @param b_x x-component of b (is not changed)
     * @param b_y y-component of b (is not changed)
     * @param b_z z-component of b (is not changed)
     */
    public static void add_b_to_a(
        float[] a_x,
        float[] a_y,
        float[] a_z,
        float[] b_x,
        float[] b_y,
        float[] b_z
    ) {
        for (int i = 0; i < a_x.length; i++) {
            a_x[i] += b_x[i];
            a_y[i] += b_y[i];
            a_z[i] += b_z[i];
        }
    }

    /**
     * Fills v for fixed molecules
     * (No checks are performed)
     * 
     * @param aCurrentTimeStep Current time step
     * @param aMoleculeVelocityFixationInfos Molecule velocity fixation infos
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     */
    public static void fill_v_forFixedMolecules(
        int aCurrentTimeStep,
        MoleculeVelocityFixationInfo[] aMoleculeVelocityFixationInfos,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z
    ) {
        for (int i = 0; i < aMoleculeVelocityFixationInfos.length; i++) {
            if (aCurrentTimeStep <= aMoleculeVelocityFixationInfos[i].getMaxTimeStep()) {
                // <editor-fold defaultstate="collapsed" desc="Fill v for fixed molecules">
                if (aMoleculeVelocityFixationInfos[i].isFixedX()) {
                    Arrays.fill(
                        aV_x,
                        aMoleculeVelocityFixationInfos[i].getFirstIndex(),
                        aMoleculeVelocityFixationInfos[i].getExclusiveLastIndex(),
                        aMoleculeVelocityFixationInfos[i].getVelocityX()
                    );
                }
                if (aMoleculeVelocityFixationInfos[i].isFixedY()) {
                    Arrays.fill(
                        aV_y,
                        aMoleculeVelocityFixationInfos[i].getFirstIndex(),
                        aMoleculeVelocityFixationInfos[i].getExclusiveLastIndex(),
                        aMoleculeVelocityFixationInfos[i].getVelocityY()
                    );
                }
                if (aMoleculeVelocityFixationInfos[i].isFixedZ()) {
                    Arrays.fill(
                        aV_z,
                        aMoleculeVelocityFixationInfos[i].getFirstIndex(),
                        aMoleculeVelocityFixationInfos[i].getExclusiveLastIndex(),
                        aMoleculeVelocityFixationInfos[i].getVelocityZ()
                    );
                }
                // </editor-fold>
            }
        }
    }
    
    /**
     * Copies rOld to r for fixed molecules
     * (No checks are performed)
     * 
     * @param aCurrentTimeStep Current time step
     * @param aMoleculeFixationInfos Molecule fixation infos
     * @param aR_x Current x-position of particle in simulation box (may be changed)
     * @param aR_y Current y-position of particle in simulation box (may be changed)
     * @param aR_z Current z-position of particle in simulation box (may be changed)
     * @param aRold_x Old x-position of particle in simulation box (is NOT changed)
     * @param aRold_y Old y-position of particle in simulation box (is NOT changed)
     * @param aRold_z Old z-position of particle in simulation box (is NOT changed)
     */
    public static void copy_rOld_to_r_forFixedMolecules(
        int aCurrentTimeStep,
        MoleculeFixationInfo[] aMoleculeFixationInfos,
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aRold_x,
        float[] aRold_y,
        float[] aRold_z
    ) {
        for (int i = 0; i < aMoleculeFixationInfos.length; i++) {
            if (aCurrentTimeStep <= aMoleculeFixationInfos[i].getMaxTimeStep()) {
                // <editor-fold defaultstate="collapsed" desc="Copy rOld to r for fixed molecules">
                if (aMoleculeFixationInfos[i].isFixedX()) {
                    System.arraycopy(
                        aRold_x, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aR_x, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aMoleculeFixationInfos[i].getLength()
                    );
                }
                if (aMoleculeFixationInfos[i].isFixedY()) {
                    System.arraycopy(
                        aRold_y, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aR_y, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aMoleculeFixationInfos[i].getLength()
                    );
                }
                if (aMoleculeFixationInfos[i].isFixedZ()) {
                    System.arraycopy(
                        aRold_z, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aR_z, 
                        aMoleculeFixationInfos[i].getFirstIndex(), 
                        aMoleculeFixationInfos[i].getLength()
                    );
                }
                // </editor-fold>
            }
        }
    }
    
    /**
     * Corrects r and v for molecule boundaries
     * (No checks are performed)
     * 
     * @param aCurrentTimeStep Current time step
     * @param aMoleculeBoundaryInfos Molecule boundary infos
     * @param aR_x Current x-position of particle in simulation box (may be changed)
     * @param aR_y Current y-position of particle in simulation box (may be changed)
     * @param aR_z Current z-position of particle in simulation box (may be changed)
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     */
    public static void correct_r_and_v_forMoleculeBoundaries(
        int aCurrentTimeStep,
        MoleculeBoundaryInfo[] aMoleculeBoundaryInfos,
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z
    ) {
        for (MoleculeBoundaryInfo tmpMoleculeBoundaryInfo: aMoleculeBoundaryInfos) {
            if (aCurrentTimeStep <= tmpMoleculeBoundaryInfo.getMaxTimeStep()) {
                if (tmpMoleculeBoundaryInfo.isBoundaryX()) {
                    Utils.correct_r_and_v_ComponentforMoleculeBoundaries(
                        aR_x,
                        aV_x,
                        tmpMoleculeBoundaryInfo.getXmin(),
                        tmpMoleculeBoundaryInfo.getXmax(),
                        tmpMoleculeBoundaryInfo.getFirstIndex(),
                        tmpMoleculeBoundaryInfo.getExclusiveLastIndex()
                    );
                }
                if (tmpMoleculeBoundaryInfo.isBoundaryY()) {
                    Utils.correct_r_and_v_ComponentforMoleculeBoundaries(
                        aR_y,
                        aV_y,
                        tmpMoleculeBoundaryInfo.getYmin(),
                        tmpMoleculeBoundaryInfo.getYmax(),
                        tmpMoleculeBoundaryInfo.getFirstIndex(),
                        tmpMoleculeBoundaryInfo.getExclusiveLastIndex()
                    );
                }
                if (tmpMoleculeBoundaryInfo.isBoundaryZ()) {
                    Utils.correct_r_and_v_ComponentforMoleculeBoundaries(
                        aR_z,
                        aV_z,
                        tmpMoleculeBoundaryInfo.getZmin(),
                        tmpMoleculeBoundaryInfo.getZmax(),
                        tmpMoleculeBoundaryInfo.getFirstIndex(),
                        tmpMoleculeBoundaryInfo.getExclusiveLastIndex()
                    );
                }
            }
        }
    }
    
    /**
     * Corrects r and v for molecule spheres
     * (No checks are performed)
     * 
     * @param aCurrentTimeStep Current time step
     * @param aMoleculeSphereInfos Molecule sphere infos
     * @param aR_x Current x-position of particle in simulation box (may be changed)
     * @param aR_y Current y-position of particle in simulation box (may be changed)
     * @param aR_z Current z-position of particle in simulation box (may be changed)
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     */
    public static void correct_r_and_v_forMoleculeSpheres(
        int aCurrentTimeStep,
        MoleculeSphereInfo[] aMoleculeSphereInfos,
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z
    ) {
        for (MoleculeSphereInfo tmpMoleculeSphereInfo: aMoleculeSphereInfos) {
            if (aCurrentTimeStep <= tmpMoleculeSphereInfo.getMaxTimeStep()) {
                int tmpFirstIndex = tmpMoleculeSphereInfo.getFirstIndex();
                int tmpExclusiveLastIndex = tmpMoleculeSphereInfo.getExclusiveLastIndex();
                for (int i = tmpFirstIndex; i < tmpExclusiveLastIndex; i++) {
                    float tmpDeltaX = aR_x[i] - tmpMoleculeSphereInfo.getSphereCenterX();
                    float tmpDeltaY = aR_y[i] - tmpMoleculeSphereInfo.getSphereCenterY();
                    float tmpDeltaZ = aR_z[i] - tmpMoleculeSphereInfo.getSphereCenterZ();
                    float tmpDistanceSquare = tmpDeltaX * tmpDeltaX + tmpDeltaY * tmpDeltaY + tmpDeltaZ * tmpDeltaZ;
                    boolean tmpIsInsideSphere = tmpDistanceSquare <= tmpMoleculeSphereInfo.getSphereRadiusSquare();
                    if (tmpMoleculeSphereInfo.isExclusiveSphere() && tmpIsInsideSphere) {
                        // <editor-fold defaultstate="collapsed" desc="Exclusive sphere (no particle allowed INSIDE sphere)">
                        float tmpDistance = (float) FastMath.sqrt(tmpDistanceSquare);
                        // NOTE: Unit vector tmpR_new_0 points from the center of the sphere in radial direction towards aR[i]
                        float tmpR_new_0_x = tmpDeltaX/tmpDistance;
                        float tmpR_new_0_y = tmpDeltaY/tmpDistance;
                        float tmpR_new_0_z = tmpDeltaZ/tmpDistance;
                        // New point outside sphere in radial direction
                        float tmpNewLength = tmpMoleculeSphereInfo.getSphereDiameter() - tmpDistance;
                        aR_x[i] = tmpMoleculeSphereInfo.getSphereCenterX() + tmpR_new_0_x * tmpNewLength;
                        aR_y[i] = tmpMoleculeSphereInfo.getSphereCenterY() + tmpR_new_0_y * tmpNewLength;
                        aR_z[i] = tmpMoleculeSphereInfo.getSphereCenterZ() + tmpR_new_0_z * tmpNewLength;
                        // New velocity vector in radial direction outside sphere with old length value
                        // float tmpV = (float) FastMath.sqrt(aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i]);
                        // aV_x[i] = tmpR_new_0_x * tmpV;
                        // aV_y[i] = tmpR_new_0_y * tmpV;
                        // aV_z[i] = tmpR_new_0_z * tmpV;
                        // Reverse velocity vector
                        aV_x[i] = -aV_x[i];
                        aV_y[i] = -aV_y[i];
                        aV_z[i] = -aV_z[i];
                        // </editor-fold>
                    } else if (!tmpMoleculeSphereInfo.isExclusiveSphere() && !tmpIsInsideSphere) {
                        // <editor-fold defaultstate="collapsed" desc="Inclusive sphere (no particle allowed OUTSIDE sphere)">
                        float tmpDistance = (float) FastMath.sqrt(tmpDistanceSquare);
                        // NOTE: Unit vector tmpR_new_0 points from the center of the sphere in radial direction towards aR[i]
                        float tmpR_new_0_x = tmpDeltaX/tmpDistance;
                        float tmpR_new_0_y = tmpDeltaY/tmpDistance;
                        float tmpR_new_0_z = tmpDeltaZ/tmpDistance;
                        // New point inside sphere in radial direction
                        float tmpNewLength = tmpMoleculeSphereInfo.getSphereDiameter() - tmpDistance;
                        aR_x[i] = tmpMoleculeSphereInfo.getSphereCenterX() + tmpR_new_0_x * tmpNewLength;
                        aR_y[i] = tmpMoleculeSphereInfo.getSphereCenterY() + tmpR_new_0_y * tmpNewLength;
                        aR_z[i] = tmpMoleculeSphereInfo.getSphereCenterZ() + tmpR_new_0_z * tmpNewLength;
                        // New velocity vector in negative radial direction inside sphere with old length value
                        // float tmpV = (float) FastMath.sqrt(aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i]);
                        // aV_x[i] = -tmpR_new_0_x * tmpV;
                        // aV_y[i] = -tmpR_new_0_y * tmpV;
                        // aV_z[i] = -tmpR_new_0_z * tmpV;
                        // Reverse velocity vector
                        aV_x[i] = -aV_x[i];
                        aV_y[i] = -aV_y[i];
                        aV_z[i] = -aV_z[i];
                        // </editor-fold>
                    }
                }
            }
        }
    }
    
    /**
     * Calculates v with f
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aF_x Current x-components of particle forces
     * @param aF_y Current y-components of particle forces
     * @param aF_z Current z-components of particle forces
     * @param aDpdMasses DPD masses of particles
     * @param aTimeStepLengthHalf Half of time step length
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    public static void calculate_v_with_f(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aF_x,
        float[] aF_y,
        float[] aF_z,
        float[] aDpdMasses,
        float aTimeStepLengthHalf,
        boolean anIsDpdUnitMass
    ) {
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] += aTimeStepLengthHalf * aF_x[i];
                aV_y[i] += aTimeStepLengthHalf * aF_y[i];
                aV_z[i] += aTimeStepLengthHalf * aF_z[i];
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] += aTimeStepLengthHalf * aF_x[i] / aDpdMasses[i];
                aV_y[i] += aTimeStepLengthHalf * aF_y[i] / aDpdMasses[i];
                aV_z[i] += aTimeStepLengthHalf * aF_z[i] / aDpdMasses[i];
            }
        }
    }

    /**
     * Adds gravitational force
     * (No checks are performed)
     * 
     * @param aGravitationalAcceleration Gravitational acceleration
     * @param aF_x Current x-components of particle forces
     * @param aF_y Current y-components of particle forces
     * @param aF_z Current z-components of particle forces
     * @param aDpdMasses DPD masses of particles
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    public static void addGravitationalForceTo_f(
        GravitationalAcceleration aGravitationalAcceleration,
        float[] aF_x,
        float[] aF_y,
        float[] aF_z,
        float[] aDpdMasses,
        boolean anIsDpdUnitMass        
    ) {
        if (anIsDpdUnitMass) {
            if (aGravitationalAcceleration.getGravitationalAccelerationX() != 0.0f) {
                float tmpGravitationalAccelerationX = aGravitationalAcceleration.getGravitationalAccelerationX();
                for (int i = 0; i < aF_x.length; i++) {
                    aF_x[i] += tmpGravitationalAccelerationX;
                }
            }
            if (aGravitationalAcceleration.getGravitationalAccelerationY() != 0.0f) {
                float tmpGravitationalAccelerationY = aGravitationalAcceleration.getGravitationalAccelerationY();
                for (int i = 0; i < aF_y.length; i++) {
                    aF_y[i] += tmpGravitationalAccelerationY;
                }
            }
            if (aGravitationalAcceleration.getGravitationalAccelerationZ() != 0.0f) {
                float tmpGravitationalAccelerationZ = aGravitationalAcceleration.getGravitationalAccelerationZ();
                for (int i = 0; i < aF_z.length; i++) {
                    aF_z[i] += tmpGravitationalAccelerationZ;
                }
            }
        } else {
            if (aGravitationalAcceleration.getGravitationalAccelerationX() != 0.0f) {
                float tmpGravitationalAccelerationX = aGravitationalAcceleration.getGravitationalAccelerationX();
                for (int i = 0; i < aF_x.length; i++) {
                    aF_x[i] += aDpdMasses[i] * tmpGravitationalAccelerationX;
                }
            }
            if (aGravitationalAcceleration.getGravitationalAccelerationY() != 0.0f) {
                float tmpGravitationalAccelerationY = aGravitationalAcceleration.getGravitationalAccelerationY();
                for (int i = 0; i < aF_y.length; i++) {
                    aF_y[i] += aDpdMasses[i] * tmpGravitationalAccelerationY;
                }
            }
            if (aGravitationalAcceleration.getGravitationalAccelerationZ() != 0.0f) {
                float tmpGravitationalAccelerationZ = aGravitationalAcceleration.getGravitationalAccelerationZ();
                for (int i = 0; i < aF_z.length; i++) {
                    aF_z[i] += aDpdMasses[i] * tmpGravitationalAccelerationZ;
                }
            }
        }
    }

    /**
     * Adds single molecule acceleration to f.
     * (No checks are performed)
     * 
     * @param aMoleculeAccelerationInfo Molecule acceleration info
     * @param aF_x Current x-components of particle forces
     * @param aF_y Current y-components of particle forces
     * @param aF_z Current z-components of particle forces
     * @param aDpdMasses DPD masses of particles
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    public static void addMoleculeAccelerationTo_f(
        MoleculeAccelerationInfo aMoleculeAccelerationInfo,
        float[] aF_x,
        float[] aF_y,
        float[] aF_z,
        float[] aDpdMasses,
        boolean anIsDpdUnitMass        
    ) {
        if (anIsDpdUnitMass) {
            if (aMoleculeAccelerationInfo.getAccelerationX() != 0.0f) {
                float tmpAccelerationX = aMoleculeAccelerationInfo.getAccelerationX();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_x[i] += tmpAccelerationX;
                }
            }
            if (aMoleculeAccelerationInfo.getAccelerationY() != 0.0f) {
                float tmpAccelerationY = aMoleculeAccelerationInfo.getAccelerationY();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_y[i] += tmpAccelerationY;
                }
            }
            if (aMoleculeAccelerationInfo.getAccelerationZ() != 0.0f) {
                float tmpAccelerationZ = aMoleculeAccelerationInfo.getAccelerationZ();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_z[i] += tmpAccelerationZ;
                }
            }
        } else {
            if (aMoleculeAccelerationInfo.getAccelerationX() != 0.0f) {
                float tmpAccelerationX = aMoleculeAccelerationInfo.getAccelerationX();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_x[i] += aDpdMasses[i] * tmpAccelerationX;
                }
            }
            if (aMoleculeAccelerationInfo.getAccelerationY() != 0.0f) {
                float tmpAccelerationY = aMoleculeAccelerationInfo.getAccelerationY();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_y[i] += aDpdMasses[i] * tmpAccelerationY;
                }
            }
            if (aMoleculeAccelerationInfo.getAccelerationZ() != 0.0f) {
                float tmpAccelerationZ = aMoleculeAccelerationInfo.getAccelerationZ();
                for (int i = aMoleculeAccelerationInfo.getFirstIndex(); i < aMoleculeAccelerationInfo.getExclusiveLastIndex(); i++) {
                    aF_z[i] += aDpdMasses[i] * tmpAccelerationZ;
                }
            }
        }
    }
    
    /**
     * Calculates r with v
     * (No checks are performed)
     * 
     * @param aR_x Current x-position of particle in simulation box (may be changed)
     * @param aR_y Current y-position of particle in simulation box (may be changed)
     * @param aR_z Current z-position of particle in simulation box (may be changed)
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aTimeStepLength Time step length
     */
    public static void calculate_r_with_v(
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float aTimeStepLength
    ) {
        for (int i = 0; i < aR_x.length; i++) {
            aR_x[i] += aV_x[i] * aTimeStepLength;
            aR_y[i] += aV_y[i] * aTimeStepLength;
            aR_z[i] += aV_z[i] * aTimeStepLength;
        }
    }
    
    /**
     * Calculates kinetic energy Ukin
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aDpdMasses DPD masses of particles
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     * @return Kinetic energy Ukin
     */
    public static float calculateUkin(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aDpdMasses,
        boolean anIsDpdUnitMass
    ) {
        // NOTE: tmpUkinDouble is twice the kinetic energy
        float tmpUkinDouble = 0.0f;
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                tmpUkinDouble += aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i];
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                tmpUkinDouble += aDpdMasses[i] * (aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i]);
            }
        }
        return tmpUkinDouble * A_HALF;
    }

    /**
     * Calculates temperature (i.e. kT, k: Boltzmann constant)
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities
     * @param aV_y Current y-components of particle velocities
     * @param aV_z Current z-components of particle velocities
     * @param aDpdMasses DPD masses of particles
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     * @return Temperature (i.e. kT, k: Boltzmann constant)
     */
    public static float getTemperature(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aDpdMasses,
        boolean anIsDpdUnitMass
    ) {
        return Utils.calculateUkin(aV_x, aV_y, aV_z, aDpdMasses, anIsDpdUnitMass) / (ONE_A_HALF * (float) aV_x.length);
    }
    
    /**
     * Calculates mean magnitude of vector
     * (No checks are performed)
     * 
     * @param aX x-component of vector
     * @param aY y-component of vector
     * @param aZ z-component of vector
     * @return Mean magnitude
     */
    public static float calculateMeanMagnitude(
        float[] aX,
        float[] aY,
        float[] aZ
    ) {
        float tmpMeanMagnitude = 0.0f;
        for (int i = 0; i < aX.length; i++) {
            tmpMeanMagnitude += (float) FastMath.sqrt(aX[i] * aX[i] + aY[i] * aY[i] + aZ[i] * aZ[i]);
        }
        return tmpMeanMagnitude / (float) aX.length;
    }
    
    /**
     * Calculates min, mean and max magnitude of vector components if magnitude
     * is not zero
     * (No checks are performed)
     * 
     * @param aX x-component of vector
     * @param aY y-component of vector
     * @param aZ z-component of vector
     * @return Array of length 3 with min (index 0), mean (index 1) and 
     * max (index 2) magnitude of (non-zero magnitude) vector components
     */
    public static float[] calculateMinMeanMax(
        float[] aX,
        float[] aY,
        float[] aZ
    ) {
        float tmpSum = 0.0f;
        float tmpMin = Float.MAX_VALUE;
        float tmpMax = 0.0f;
        int tmpCounter = 0;
        for (int i = 0; i < aX.length; i++) {
            if (aX[i] != 0.0f || aY[i] != 0.0f || aZ[i] != 0.0f) {
                tmpCounter++;
                float tmpMagnitude = (float) FastMath.sqrt(aX[i] * aX[i] + aY[i] * aY[i] + aZ[i] * aZ[i]);
                tmpSum += tmpMagnitude;
                if (tmpMagnitude < tmpMin) {
                    tmpMin = tmpMagnitude;
                }
                if (tmpMagnitude > tmpMax) {
                    tmpMax = tmpMagnitude;
                }
            }
        }
        if (tmpSum == 0.0f) {
            return null;
        } else {
            return new float[]{tmpMin, tmpSum / (float) tmpCounter, tmpMax};
        }
    }

    /**
     * Scales v
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aDpdMasses DPD masses of particles
     * @param aTemperature Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     * @return Velocity scale factor
     */
    public static float scale_v(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aDpdMasses,
        float aTemperature,
        boolean anIsDpdUnitMass
    ) {
        // Remove excess momentum first
        Utils.removeExcessMomentum(
            aV_x, 
            aV_y, 
            aV_z, 
            aDpdMasses, 
            anIsDpdUnitMass
        );
        
        float tmpVelocityScaleFactor = 
            Utils.getVelocityScaleFactor(
                aV_x, 
                aV_y, 
                aV_z, 
                aDpdMasses, 
                aTemperature, 
                anIsDpdUnitMass
            );
        for (int i = 0; i < aV_x.length; i++) {
            aV_x[i] *= tmpVelocityScaleFactor;
            aV_y[i] *= tmpVelocityScaleFactor;
            aV_z[i] *= tmpVelocityScaleFactor;
        }
        return tmpVelocityScaleFactor;
    }

    /**
     * Returns velocity scale factor
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aDpdMasses DPD masses of particles
     * @param aTemperature Temperature in DPD units (i.e. kT fractions, k: Boltzmann constant)
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     * @return Velocity scale factor
     */
    public static float getVelocityScaleFactor(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aDpdMasses,
        float aTemperature,
        boolean anIsDpdUnitMass
    ) {
        // NOTE: tmpUkinDouble is twice the real kinetic energy
        float tmpUkinDouble = 0.0f;
        if (anIsDpdUnitMass) {
            for (int i = 0; i < aV_x.length; i++) {
                tmpUkinDouble += aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i];
            }
        } else {
            for (int i = 0; i < aV_x.length; i++) {
                tmpUkinDouble += aDpdMasses[i] * (aV_x[i] * aV_x[i] + aV_y[i] * aV_y[i] + aV_z[i] * aV_z[i]);
            }
        }
        float tmpDegreesOfFreedom = THREE * (float) aV_x.length - THREE;
        return (float) FastMath.sqrt(tmpDegreesOfFreedom * aTemperature/tmpUkinDouble);
    }
    
    /**
     * Removes excess momentum from velocities
     * (No checks are performed)
     * 
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aDpdMasses DPD masses of particles
     * @param anIsDpdUnitMass Flag for use of DPD unit masses. True : DPD masses of all particles are set to 1, False: ...
     */
    public static void removeExcessMomentum(
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        float[] aDpdMasses,
        boolean anIsDpdUnitMass
    ) {
        if (anIsDpdUnitMass) {
            float tmpExcessMomentum_x = 0.0f;
            float tmpExcessMomentum_y = 0.0f;
            float tmpExcessMomentum_z = 0.0f;
            for (int i = 0; i < aV_x.length; i++) {
                tmpExcessMomentum_x += aV_x[i];
                tmpExcessMomentum_y += aV_y[i];
                tmpExcessMomentum_z += aV_z[i];
            }
            float tmpExcessMomentumPerParticle_x = tmpExcessMomentum_x / (float) aV_x.length;
            float tmpExcessMomentumPerParticle_y = tmpExcessMomentum_y / (float) aV_x.length;
            float tmpExcessMomentumPerParticle_z = tmpExcessMomentum_z / (float) aV_x.length;
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] -= tmpExcessMomentumPerParticle_x;
                aV_y[i] -= tmpExcessMomentumPerParticle_y;
                aV_z[i] -= tmpExcessMomentumPerParticle_z;
            }
        } else {
            float tmpExcessMomentum_x = 0.0f;
            float tmpExcessMomentum_y = 0.0f;
            float tmpExcessMomentum_z = 0.0f;
            for (int i = 0; i < aV_x.length; i++) {
                tmpExcessMomentum_x += aDpdMasses[i] * aV_x[i];
                tmpExcessMomentum_y += aDpdMasses[i] * aV_y[i];
                tmpExcessMomentum_z += aDpdMasses[i] * aV_z[i];
            }
            float tmpExcessMomentumPerParticle_x = tmpExcessMomentum_x / (float) aV_x.length;
            float tmpExcessMomentumPerParticle_y = tmpExcessMomentum_y / (float) aV_x.length;
            float tmpExcessMomentumPerParticle_z = tmpExcessMomentum_z / (float) aV_x.length;
            for (int i = 0; i < aV_x.length; i++) {
                aV_x[i] -= tmpExcessMomentumPerParticle_x / aDpdMasses[i];
                aV_y[i] -= tmpExcessMomentumPerParticle_y / aDpdMasses[i];
                aV_z[i] -= tmpExcessMomentumPerParticle_z / aDpdMasses[i];
            }
        }
    }
    
    /**
     * Corrects position difference for threshold and periodic boundary 
     * 
     * @param aPositionDifference Difference between two positions
     * @param anIsPeriodicBoundaryAlongAxis True: Periodic boundary along axis, false: Otherwise
     * @param aBoxAxisHalfLength Half length of box along axis
     * @param aBoxAxisLength Length of box along axis
     * @param aBoxAxisNegativeHalfLength Negative half length of box along axis
     * @return Correct position difference
     */
    public static float correctPositionDifference(
        float aPositionDifference,
        boolean anIsPeriodicBoundaryAlongAxis,
        float aBoxAxisHalfLength,
        float aBoxAxisLength,
        float aBoxAxisNegativeHalfLength
    ) {
        if (FastMath.abs(aPositionDifference) < Constants.TINY_THRESHOLD) {
            if (aPositionDifference < 0.0f) {
                return -Constants.TINY_THRESHOLD;
            } else {
                return Constants.TINY_THRESHOLD;
            }
        } else {
            if (anIsPeriodicBoundaryAlongAxis) {
                if (aPositionDifference > aBoxAxisHalfLength) {
                    aPositionDifference -= aBoxAxisLength;
                } else if (aPositionDifference < aBoxAxisNegativeHalfLength) {
                    aPositionDifference += aBoxAxisLength;
                }
            }
            return aPositionDifference;
        }
    }

    /**
     * Corrects r and v
     * (No checks are performed)
     * 
     * @param aR_x x-position of particle in simulation box (may be changed)
     * @param aR_y y-position of particle in simulation box (may be changed)
     * @param aR_z z-position of particle in simulation box (may be changed)
     * @param aV_x Current x-components of particle velocities (may be changed)
     * @param aV_y Current y-components of particle velocities (may be changed)
     * @param aV_z Current z-components of particle velocities (may be changed)
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries 
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     */
    public static void correct_r_and_v(
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        float[] aV_x,
        float[] aV_y,
        float[] aV_z,
        BoxSize aBoxSize,
        PeriodicBoundaries aPeriodicBoundaries,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        Utils.correct_r_and_v_Component(
            aR_x,
            aV_x,
            aPeriodicBoundaries.isAlongX(),
            aBoxSize.getXLength(),
            aBoxSize.getXDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
        Utils.correct_r_and_v_Component(
            aR_y,
            aV_y,
            aPeriodicBoundaries.isAlongY(),
            aBoxSize.getYLength(),
            aBoxSize.getYDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
        Utils.correct_r_and_v_Component(
            aR_z,
            aV_z,
            aPeriodicBoundaries.isAlongZ(),
            aBoxSize.getZLength(),
            aBoxSize.getZDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
    };

    /**
     * Calculates initial potential energy minimization steps
     * 
     * @param aSimulationLogger Simulation logger
     * @param aParameters Parameters
     * @param aFactory Factory
     * @param aConservativeForceAccumulator Conservative force accumulator (is shutdown after operation)
     * @param aPotentialAccumulator Potential accumulator
     * @param aSimulationOutput Simulation output
     * @param aParticlePositionPool Particle position pool
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     * @throws IllegalArgumentException Thrown if argument is illegal
     */
    public static void calculateInitialPotentialEnergyMinimizationSteps(
        ILogger aSimulationLogger,
        Parameters aParameters,
        Factory aFactory,
        ForceAccumulator aConservativeForceAccumulator,
        PotentialAccumulator aPotentialAccumulator,
        IOutput aSimulationOutput,
        ParticlePositionPool aParticlePositionPool,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        // <editor-fold defaultstate="collapsed" desc="Checks">
        if (aSimulationLogger == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aSimulationLogger is null.");
        }
        if (aParameters == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aParameters is null.");
        }
        if (aFactory == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aFactory is null.");
        }
        if (aConservativeForceAccumulator == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aConservativeForceAccumulator is null.");
        }
        if (aPotentialAccumulator == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aPotentialAccumulator is null.");
        }
        if (aSimulationOutput == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aSimulationOutput is null.");
        }
        if (aParticlePositionPool == null) {
            throw new IllegalArgumentException("Utils.calculateInitialPotentialEnergyMinimizationSteps: aParticlePositionPool is null.");
        }
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        long tmpId = aSimulationLogger.getId();
        aSimulationLogger.appendMethodCallStart("Utils.calculateInitialPotentialEnergyMinimizationSteps", tmpId);
        // </editor-fold>
        ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        SimulationDescription tmpSimulationDescription = aParameters.getSimulationDescription();
        ChemicalSystemDescription tmpChemicalSystemDescription = aParameters.getChemicalSystemDescription();
        float[] tmpR_x = tmpParticleArrays.getR_x();
        float[] tmpR_y = tmpParticleArrays.getR_y();
        float[] tmpR_z = tmpParticleArrays.getR_z();
        float[] tmpRold_x = tmpParticleArrays.getRold_x();
        float[] tmpRold_y = tmpParticleArrays.getRold_y();
        float[] tmpRold_z = tmpParticleArrays.getRold_z();
        float[] tmpF_x = tmpParticleArrays.getF_x();
        float[] tmpF_y = tmpParticleArrays.getF_y();
        float[] tmpF_z = tmpParticleArrays.getF_z();
        BoxSize tmpBoxSize = tmpChemicalSystemDescription.getBoxSize();
        PeriodicBoundaries tmpPeriodicBoundaries = tmpSimulationDescription.getPeriodicBoundaries();
        // Step length in DPD units derived from the minimum box side length
        float tmpStepLength = 0.00001f * aConservativeForceAccumulator.getParticlePairDpdForceCalculator().getBoxSize().getMinimumLength();
        aPotentialAccumulator.accumulatePotentials(
            tmpParticleArrays.getBondChunkArraysList(),
            tmpR_x,
            tmpR_y,
            tmpR_z,
            aParameters
        );
        float tmpUpotMin =aPotentialAccumulator.getExtendedAdderGroup().getPotentialEnergyAdder().getSum();
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Upot(START)       = " + String.valueOf(tmpUpotMin));
        // </editor-fold>
        for (int i = 0; i < tmpSimulationDescription.getInitialPotentialEnergyMinimizationStepNumber(); i++) {
            // Save x,y,z-positions
            for (int k = 0; k < tmpR_x.length; k++) {
                tmpRold_x[k] = tmpR_x[k];
                tmpRold_y[k] = tmpR_y[k];
                tmpRold_z[k] = tmpR_z[k];
            }
            aConservativeForceAccumulator.accumulate_f_and_fTwo(
                tmpParticleArrays.getBondChunkArraysList(),
                tmpR_x,
                tmpR_y,
                tmpR_z,
                aParameters
            );
            float tmpMaximum = 0.0f;
            for (int k = 0; k < tmpR_x.length; k++) {
                tmpMaximum = (float) FastMath.max((float) FastMath.abs(tmpF_x[k]), (float) FastMath.max((float) FastMath.abs(tmpF_y[k]), (float) FastMath.abs(tmpF_z[k])));
            }
            float tmpFactor = tmpStepLength / tmpMaximum;
            for (int k = 0; k < tmpR_x.length; k++) {
                // Displace particles
                tmpR_x[k] += tmpFactor * tmpF_x[k];
                tmpR_y[k] += tmpFactor * tmpF_y[k];
                tmpR_z[k] += tmpFactor * tmpF_z[k];
            }
            Utils.correct_r(
                tmpR_x,
                tmpR_y,
                tmpR_z,
                tmpBoxSize,
                tmpPeriodicBoundaries,
                aMaximumNumberOfPositionCorrectionTrials
            );
            
            aPotentialAccumulator.accumulatePotentials(
                tmpParticleArrays.getBondChunkArraysList(),
                tmpR_x,
                tmpR_y,
                tmpR_z,
                aParameters
            );
            float tmpUpotCurrent = aPotentialAccumulator.getExtendedAdderGroup().getPotentialEnergyAdder().getSum();
            
            // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
            aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Minimization step = " + String.valueOf(i + 1) + ", Timestep Length = " + String.valueOf(tmpStepLength));
            // </editor-fold>
            if (tmpUpotCurrent < tmpUpotMin) {
                tmpStepLength *= 2.0f;
                tmpUpotMin = tmpUpotCurrent;
                // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Minimization step = " + String.valueOf(i + 1) + ", Upot(Min)       = " + String.valueOf(tmpUpotMin));
                // </editor-fold>
            } else {
                // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Minimization step = " + String.valueOf(i + 1) + ", Upot(Current)   = " + String.valueOf(tmpUpotCurrent));
                // </editor-fold>
                // Restore x,y,z-positions
                for (int k = 0; k < tmpR_x.length; k++) {
                    tmpR_x[k] = tmpRold_x[k];
                    tmpR_y[k] = tmpRold_y[k];
                    tmpR_z[k] = tmpRold_z[k];
                }
                if (tmpStepLength > 1E-12) {
                    tmpStepLength *= 0.5;
                } else {
                    // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
                    aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Minimization step = " + String.valueOf(i + 1) + ", Break");
                    // </editor-fold>
                    break;
                }
            }

            if (tmpSimulationDescription.isInitialPotentialEnergyMinimizationStepOutput()) {
                aSimulationOutput.setMinimizationStepParticlePositions(i + 1, Utils.getParticlePositions(aParameters, aParticlePositionPool));
            }
        }
        // <editor-fold defaultstate="collapsed" desc="Intermediate results logging">
        aSimulationLogger.appendIntermediateResults("Utils.calculateInitialPotentialEnergyMinimizationSteps, Upot(Minimum) = " + String.valueOf(tmpUpotMin));
        // </editor-fold>
        // <editor-fold defaultstate="collapsed" desc="Method call logging">
        aSimulationLogger.appendMethodCallEnd("Utils.calculateInitialPotentialEnergyMinimizationSteps", tmpId);
        // </editor-fold>
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Particle pair related methods">
    /**
     * Returns particle token pair key.
     * (No checks are performed)
     *
     * @param aParticleToken1 Particle token 1 of pair
     * @param aParticleToken2 Particle token 2 of pair
     * @return Particle token pair key
     */
    public static String getParticleTokenPairKey(String aParticleToken1, String aParticleToken2) {
        StringBuilder tmpBuffer = new StringBuilder(aParticleToken1.length() + GENERAL_SEPARATOR.length() + aParticleToken2.length());
        // Sort particle names so that particle pair is always unique (avoids different entries for A-B and B-A)
        if (aParticleToken1.compareTo(aParticleToken2) < 0) {
            tmpBuffer.append(aParticleToken1);
            tmpBuffer.append(GENERAL_SEPARATOR);
            tmpBuffer.append(aParticleToken2);
        } else {
            tmpBuffer.append(aParticleToken2);
            tmpBuffer.append(GENERAL_SEPARATOR);
            tmpBuffer.append(aParticleToken1);
        }
        return tmpBuffer.toString();
    }

    /**
     * Returns particle index pair key.
     * (No checks are performed)
     *
     * @param aParticleIndex1 Particle index 1 of pair
     * @param aParticleIndex2 Particle index 2 of pair
     * @return Particle index pair key
     */
    public static String getParticleIndexPairKey(int aParticleIndex1, int aParticleIndex2) {
        return Utils.getParticleTokenPairKey(String.valueOf(aParticleIndex1), String.valueOf(aParticleIndex2));
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Particle position related methods">
    /**
     * Returns particle positions.
     * (No checks are performed)
     * NOTE: Particle position instances are obtained from specified particle 
     * position pool.
     * 
     * @param aParameters Parameters
     * @param aParticlePositionPool Particle position pool
     * @return Particle positions
     */
    public static ParticlePosition[] getParticlePositions(Parameters aParameters, ParticlePositionPool aParticlePositionPool) {
        ParticleArrays tmpParticleArrays = aParameters.getParticleArrays();
        String[] tmpTokens = tmpParticleArrays.getParticleTokens();
        int[] tmpMoleculeTypeIndices = tmpParticleArrays.getMoleculeTypeIndices();
        float[] tmpR_x = tmpParticleArrays.getR_x();
        float[] tmpR_y = tmpParticleArrays.getR_y();
        float[] tmpR_z = tmpParticleArrays.getR_z();
        int[] tmpMoleculeIndices = tmpParticleArrays.getMoleculeIndices();
        MoleculeTypes tmpMoleculeTypes = aParameters.getChemicalSystemDescription().getMoleculeTypes();
        ParticlePosition[] tmpParticlePositions = new ParticlePosition[aParameters.getSimulationCounts().getParticleNumber()];
        for (int i = 0; i < aParameters.getSimulationCounts().getParticleNumber(); i++) {
            ParticlePosition tmpParticlePosition = aParticlePositionPool.getParticlePosition();
            tmpParticlePosition.setPosition(
                tmpTokens[i],
                tmpMoleculeTypes.getMoleculeName(tmpMoleculeTypeIndices[i]),
                tmpR_x[i],
                tmpR_y[i],
                tmpR_z[i],
                i,
                tmpMoleculeIndices[i]
            );
            tmpParticlePositions[i] = tmpParticlePosition;
        }
        return tmpParticlePositions;
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Molecule-particles related methods">
    /**
     * Returns molecule-particle string.
     * (No checks are performed)
     * 
     * @param aParticleToken Particle token
     * @param aMoleculeName Name of molecule of particle
     * @return Molecule-particle string
     */
    public static String getMoleculeParticle(String aParticleToken, String aMoleculeName) {
        return aMoleculeName + Utils.MOLECULE_NAME_PARTICLE_TOKEN_SEPARATOR + aParticleToken;
    }
    // </editor-fold>
    // <editor-fold defaultstate="collapsed" desc="Electrostatics related methods">
    /**
     * Returns DPD electrostatics conservative force
     * (No checks are performed)
     * 
     * @param tmpParticleCharge1 Charge of particle 1
     * @param tmpParticleCharge2 Charge of particle 2
     * @param tmpElectrostatics Electrostatics object
     * @param tmpRij Rij
     * @return DPD electrostatics conservative force
     */
    public static float getElectrostaticsDpdForce(
        float tmpParticleCharge1,
        float tmpParticleCharge2,
        Electrostatics tmpElectrostatics,
        float tmpRij
    ) {
        // Coulomb term
        float tmpTermCoulomb;
        float tmpDerivativeTermCoulomb;
        if (tmpElectrostatics.getEffectiveExponent() == 2.0) {
            tmpTermCoulomb = tmpElectrostatics.getElectrostaticsCoupling() * tmpParticleCharge1 * tmpParticleCharge2 / tmpRij;
            tmpDerivativeTermCoulomb = -tmpTermCoulomb / tmpRij;
        } else {
            tmpTermCoulomb = tmpElectrostatics.getElectrostaticsCoupling() * tmpParticleCharge1 * tmpParticleCharge2 / (float) FastMath.pow(tmpRij, tmpElectrostatics.getEffectiveExponent() - 1);
            tmpDerivativeTermCoulomb = -(tmpElectrostatics.getEffectiveExponent() - 1) * tmpTermCoulomb / tmpRij;
        }
        // Charge distribution
        float tmpTermChargeDistribution = 1.0f;
        float tmpDerivativeTermChargeDistribution = 0.0f;
        if (tmpElectrostatics.getChargeDistributionType() == Electrostatics.ChargeDistributionType.ALEJANDRE) {
            float tmpDecayLengthAlejandre = tmpElectrostatics.getDecayLengthAlejandre();
            float tmpRij_div_halfLambda = TWO * tmpRij / tmpDecayLengthAlejandre;
            float tmpExp = (float) FastMath.exp(-tmpRij_div_halfLambda);
            tmpTermChargeDistribution = ONE - (ONE + A_HALF * tmpRij_div_halfLambda) * tmpExp;
            tmpDerivativeTermChargeDistribution = (ONE + tmpRij_div_halfLambda) * tmpExp / tmpDecayLengthAlejandre;
        }
        // Splitting type
        float tmpTermSplitting = 1.0f;
        float tmpDerivativeTermSplitting = 0.0f;
        if (tmpElectrostatics.getSplittingType()== Electrostatics.SplittingType.FANOURGAKIS) {
            float tmpRcutoff = tmpElectrostatics.getCutOffLength();
            float tmpRij_div_Rcutoff = tmpRij / tmpRcutoff;
            float tmpRij_div_Rcutoff_2 = tmpRij_div_Rcutoff * tmpRij_div_Rcutoff;
            float tmpRij_div_Rcutoff_4 = tmpRij_div_Rcutoff_2 * tmpRij_div_Rcutoff_2;
            float tmpRij_div_Rcutoff_5 = tmpRij_div_Rcutoff_4 * tmpRij_div_Rcutoff;
            float tmpRij_div_Rcutoff_6 = tmpRij_div_Rcutoff_4 * tmpRij_div_Rcutoff_2;
            float tmpRij_div_Rcutoff_7 = tmpRij_div_Rcutoff_6 * tmpRij_div_Rcutoff;
            tmpTermSplitting = ONE - C7_4 * tmpRij_div_Rcutoff + C21_4 * tmpRij_div_Rcutoff_5 - SEVEN * tmpRij_div_Rcutoff_6 + C5_2 * tmpRij_div_Rcutoff_7;
            tmpDerivativeTermSplitting = (C35_2 * tmpRij_div_Rcutoff_6 - FOURTY_TWO * tmpRij_div_Rcutoff_5 + C105_4 * tmpRij_div_Rcutoff_4 - C7_4) / tmpRcutoff;
        }
        // Force = - dEpot/dRij
        float tmpForce =
            -(
                tmpDerivativeTermCoulomb * tmpTermChargeDistribution * tmpTermSplitting +
                tmpTermCoulomb * tmpDerivativeTermChargeDistribution * tmpTermSplitting +
                tmpTermCoulomb * tmpTermChargeDistribution * tmpDerivativeTermSplitting
            );
        return tmpForce;
    }

    /**
     * Returns DPD electrostatics potential energy
     * (No checks are performed)
     * 
     * @param tmpParticleCharge1 Charge of particle 1
     * @param tmpParticleCharge2 Charge of particle 2
     * @param tmpElectrostatics Electrostatics object
     * @param tmpRij Rij
     * @return DPD electrostatics potential
     */
    public static float getElectrostaticsDpdPotentialEnergy(
        float tmpParticleCharge1,
        float tmpParticleCharge2,
        Electrostatics tmpElectrostatics,
        float tmpRij
    ) {
        // Coulomb term
        float tmpTermCoulomb;
        if (tmpElectrostatics.getEffectiveExponent() == 2.0) {
            tmpTermCoulomb = tmpElectrostatics.getElectrostaticsCoupling() * tmpParticleCharge1 * tmpParticleCharge2 / tmpRij;
        } else {
            tmpTermCoulomb = tmpElectrostatics.getElectrostaticsCoupling() * tmpParticleCharge1 * tmpParticleCharge2 / (float) FastMath.pow(tmpRij, tmpElectrostatics.getEffectiveExponent() - 1);
        }
        // Charge distribution
        float tmpTermChargeDistribution = 1.0f;
        if (tmpElectrostatics.getChargeDistributionType() == Electrostatics.ChargeDistributionType.ALEJANDRE) {
            float tmpDecayLengthAlejandre = tmpElectrostatics.getDecayLengthAlejandre();
            float tmpRij_div_halfLambda = TWO * tmpRij / tmpDecayLengthAlejandre;
            float tmpExp = (float) FastMath.exp(-tmpRij_div_halfLambda);
            tmpTermChargeDistribution = ONE - (ONE + A_HALF * tmpRij_div_halfLambda) * tmpExp;
        }
        // Splitting type
        float tmpTermSplitting = 1.0f;
        if (tmpElectrostatics.getSplittingType()== Electrostatics.SplittingType.FANOURGAKIS) {
            float tmpRcutoff = tmpElectrostatics.getCutOffLength();
            float tmpRij_div_Rcutoff = tmpRij / tmpRcutoff;
            float tmpRij_div_Rcutoff_2 = tmpRij_div_Rcutoff * tmpRij_div_Rcutoff;
            float tmpRij_div_Rcutoff_4 = tmpRij_div_Rcutoff_2 * tmpRij_div_Rcutoff_2;
            float tmpRij_div_Rcutoff_5 = tmpRij_div_Rcutoff_4 * tmpRij_div_Rcutoff;
            float tmpRij_div_Rcutoff_6 = tmpRij_div_Rcutoff_4 * tmpRij_div_Rcutoff_2;
            float tmpRij_div_Rcutoff_7 = tmpRij_div_Rcutoff_6 * tmpRij_div_Rcutoff;
            tmpTermSplitting = ONE - C7_4 * tmpRij_div_Rcutoff + C21_4 * tmpRij_div_Rcutoff_5 - SEVEN * tmpRij_div_Rcutoff_6 + C5_2 * tmpRij_div_Rcutoff_7;
        }
        // Epot
        float tmpPotential = tmpTermCoulomb * tmpTermChargeDistribution * tmpTermSplitting;
        return tmpPotential;
    }
    // </editor-fold>
    // </editor-fold>
    //
    // <editor-fold defaultstate="collapsed" desc="Private methods">
    /**
     * Corrects r and v component for molecule boundaries
     * (No checks are performed)
     * 
     * @param aR Current component of particle position in simulation box (may be changed)
     * @param aV Current component of particle velocities (may be changed)
     * @param aMin Minimum value
     * @param aMax Maximum value
     * @param aFirstIndex First index
     * @param anExclusiveLastIndex Exclusive last index
     */
    private static void correct_r_and_v_ComponentforMoleculeBoundaries(
        float[] aR,
        float[] aV,
        float aMin,
        float aMax,
        int aFirstIndex,
        int anExclusiveLastIndex
    ) {
        boolean tmpIsMinSmallerMax = aMin < aMax;
        boolean tmpIsMinGreaterMax = aMin > aMax;
        float tmpDoubleMin = aMin + aMin;
        float tmpDoubleMax = aMax + aMax;
        for (int i = aFirstIndex; i < anExclusiveLastIndex; i++) {
            if (tmpIsMinSmallerMax) {
                if (aR[i] < aMin) {
                    aR[i] = tmpDoubleMin - aR[i];
                    aV[i] = -aV[i];
                } else if (aR[i] > aMax) {
                    aR[i] = tmpDoubleMax - aR[i];
                    aV[i] = -aV[i];
                }
            } else if (tmpIsMinGreaterMax) {
                if (aR[i] < aMin && aR[i] > aMax) {
                    if (aMin - aR[i] < aR[i] - aMax) {
                        aR[i] = tmpDoubleMin - aR[i];
                        aV[i] = -aV[i];
                    } else {
                        aR[i] = tmpDoubleMax - aR[i];
                        aV[i] = -aV[i];
                    }
                }
            }
        }
    }

    /**
     * Corrects r
     * (No checks are performed)
     * 
     * @param aR_x x-position of particle in simulation box (may be changed)
     * @param aR_y y-position of particle in simulation box (may be changed)
     * @param aR_z z-position of particle in simulation box (may be changed)
     * @param aBoxSize Box size
     * @param aPeriodicBoundaries Periodic boundaries 
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     */
    private static void correct_r(
        float[] aR_x,
        float[] aR_y,
        float[] aR_z,
        BoxSize aBoxSize,
        PeriodicBoundaries aPeriodicBoundaries,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        Utils.correct_r_Component(
            aR_x,
            aPeriodicBoundaries.isAlongX(),
            aBoxSize.getXLength(),
            aBoxSize.getXDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
        Utils.correct_r_Component(
            aR_y,
            aPeriodicBoundaries.isAlongY(),
            aBoxSize.getYLength(),
            aBoxSize.getYDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
        Utils.correct_r_Component(
            aR_z,
            aPeriodicBoundaries.isAlongZ(),
            aBoxSize.getZLength(),
            aBoxSize.getZDoubleLength(),
            aMaximumNumberOfPositionCorrectionTrials
        );
    };

    /**
     * Corrects r component
     * (No checks are performed)
     * 
     * @param aR Current component of particle position in simulation box (may be changed)
     * @param anIsPeriodicBoundary True: Periodic boundary, false: Otherwise
     * @param aBoxSizeLength Box size length
     * @param aDoubleBoxSizeLength Double box size length
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     */
    private static void correct_r_Component(
        float[] aR,
        boolean anIsPeriodicBoundary,
        float aBoxSizeLength,
        float aDoubleBoxSizeLength,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        if (anIsPeriodicBoundary) {
            Utils.correct_r_ComponentForPeriodicBoundary(
                aR,
                aBoxSizeLength,
                aMaximumNumberOfPositionCorrectionTrials
            );
        } else {
            for (int i = 0; i < aR.length; i++) {
                for (int k = 0; k <= aMaximumNumberOfPositionCorrectionTrials; k++) {
                    if (aR[i] > aBoxSizeLength) {
                        aR[i] = aDoubleBoxSizeLength - aR[i];
                    } else if (aR[i] < 0.0f) {
                        aR[i] = -aR[i];
                    } else if (aR[i] == aBoxSizeLength) {
                        aR[i] = Utils.getSlightlyReducedValue(aBoxSizeLength);
                    } else {
                        break;
                    }
                    if (k == aMaximumNumberOfPositionCorrectionTrials) {
                        throw new IllegalStateException("Utils.correct_r_Component: A particle is outside the simulation box (aMaximumNumberOfPositionCorrectionTrials exceeded).");
                    }
                }
            }
        }
    };

    /**
     * Corrects r and v for periodic boundary
     * (No checks are performed)
     * 
     * @param aR Current component of particle position in simulation box (may be changed)
     * @param aV Current component of particle velocities (may be changed)
     * @param anIsPeriodicBoundary True: Periodic boundary, false: Otherwise
     * @param aBoxSizeLength Box size length
     * @param aDoubleBoxSizeLength Double box size length
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     */
    private static void correct_r_and_v_Component(
        float[] aR,
        float[] aV,
        boolean anIsPeriodicBoundary,
        float aBoxSizeLength,
        float aDoubleBoxSizeLength,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        if (anIsPeriodicBoundary) {
            Utils.correct_r_ComponentForPeriodicBoundary(
                aR,
                aBoxSizeLength,
                aMaximumNumberOfPositionCorrectionTrials
            );
        } else {
            for (int i = 0; i < aR.length; i++) {
                for (int k = 0; k <= aMaximumNumberOfPositionCorrectionTrials; k++) {
                    if (aR[i] > aBoxSizeLength) {
                        aR[i] = aDoubleBoxSizeLength - aR[i];
                        aV[i] = -aV[i];
                    } else if (aR[i] < 0.0f) {
                        aR[i] = -aR[i];
                        aV[i] = -aV[i];
                    } else if (aR[i] == aBoxSizeLength) {
                        aR[i] = Utils.getSlightlyReducedValue(aBoxSizeLength);
                        aV[i] = -aV[i];
                    } else {
                        break;
                    }
                    if (k == aMaximumNumberOfPositionCorrectionTrials) {
                        throw new IllegalStateException("Utils.correct_r_and_v_Component: A particle is outside the simulation box (aMaximumNumberOfPositionCorrectionTrials exceeded).");
                    }
                }
            }
        }
    };
    
    /**
     * Corrects r component for periodic boundary
     * (No checks are performed)
     * 
     * @param aR Current component of particle position in simulation box (may be changed)
     * @param anIndex Index
     * @param aBoxSizeLength Box size length
     * @param aMaximumNumberOfPositionCorrectionTrials Maximum number of position correction trials
     */
    private static void correct_r_ComponentForPeriodicBoundary(
        float[] aR,
        float aBoxSizeLength,
        int aMaximumNumberOfPositionCorrectionTrials
    ) {
        for (int i = 0; i < aR.length; i++) {
            for (int k = 0; k <= aMaximumNumberOfPositionCorrectionTrials; k++) {
                if (aR[i] > aBoxSizeLength) {
                    aR[i] -= aBoxSizeLength;
                } else if (aR[i] < 0.0f) {
                    aR[i] += aBoxSizeLength;
                } else if (aR[i] == aBoxSizeLength) {
                    aR[i] = 0.0f;
                } else {
                    break;
                }
                if (k == aMaximumNumberOfPositionCorrectionTrials) {
                    throw new IllegalStateException("Utils.correct_r_ComponentForPeriodicBoundary: A particle is outside the simulation box (aMaximumNumberOfPositionCorrectionTrials exceeded).");
                }
            }
        }
    };
    
    /**
     * Returns slightly reduced value
     * 
     * @param aValue Value to be slightly reduced
     * @return Slightly reduced value
     */
    private static float getSlightlyReducedValue(float aValue) {
        float tmpNewValue = aValue;
        int tmpCounter = 1;
        while (tmpNewValue == aValue) {
            tmpNewValue *= 1.0f - (float) tmpCounter * Constants.TINY_THRESHOLD;
            tmpCounter++;
        }
        return tmpNewValue;
    }
    // </editor-fold>
    
}
