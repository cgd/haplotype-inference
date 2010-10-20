/*
 * Copyright (c) 2010 The Jackson Laboratory
 * 
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.data;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.codec.net.URLCodec;
import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.haplotype.io.HiddenMarkovModelStateParser;
import org.jax.util.datastructure.SequenceUtilities;
import org.jax.util.io.CommonFlatFileFormat;
import org.jax.util.io.FileExtensionFilter;
import org.jax.util.io.FlatFileReader;
import org.jax.util.io.FlatFileWriter;
import org.jax.util.io.IllegalFormatException;

/**
 * TODO finish implementing me!
 * 
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class BinaryMultiGroupHaplotypeDataSource implements MultiGroupHaplotypeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 7948485682227419342L;
    
    private static final URLCodec URL_CODEC = new URLCodec();
    public static final FileExtensionFilter HMM_STATE_STREAM_FILTER =
        new FileExtensionFilter("hmmstates");
    private static final String CHROMOSOME_DIR_PREFIX = "chr";
    private static final String STRAIN_NAMES_FILENAME = "strain-names.csv";
    
    private final String name;

    private final File directory;
    
    /**
     * Constructor
     * @param name
     *          the name of this genome data source
     * @param directory
     *          the directory
     */
    public BinaryMultiGroupHaplotypeDataSource(
            String name,
            File directory)
    {
        this.name = name;
        this.directory = directory;
    }
    
    /**
     * {@inheritDoc}
     */
    public List<MultiPartitionedInterval> getHaplotypeData(
            int chromosome,
            Set<String> strainsToAccept)
    {
        // TODO Auto-generated method stub
        return null;
    }
    
    /**
     * Getter for the available strains
     * @return
     *          the available strains
     */
    public Set<String> getAvailableStrains()
    {
        try
        {
            return BinaryMultiGroupHaplotypeDataSource.readStrainsFromCsvData(
                    this.directory);
        }
        catch(RuntimeException ex)
        {
            throw ex;
        }
        catch(Exception ex)
        {
            // the caller should not have to guard against exceptions here
            throw new RuntimeException(ex);
        }
    }
    
    /**
     * Getter for the available chromosomes
     * @return
     *          the available chromosomes
     */
    public int[] getAvailableChromosomes()
    {
        return getAvailableBinaryDataChromosomes(this.directory);
    }
    
    /**
     * Getter for the name
     * @return the name
     */
    public String getName()
    {
        return this.name;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return this.getName();
    }
    
    /**
     * Read in the strains from the given comma-separated data
     * @param directory
     *          the directory to use
     * @return
     *          the strains
     * @throws IOException
     *          if file IO fails
     * @throws IllegalFormatException
     *          if the file format is bad
     */
    public static Set<String> readStrainsFromCsvData(
            File directory) throws IOException, IllegalFormatException
    {
        File csvFile = new File(directory, STRAIN_NAMES_FILENAME);
        FlatFileReader ffr = new FlatFileReader(
                new FileReader(csvFile),
                CommonFlatFileFormat.CSV_UNIX);
        String[] strainRow = ffr.readRow();
        
        return new HashSet<String>(Arrays.asList(strainRow));
    }
    
    /**
     * Write strains as CSV
     * @param strains
     *          the strains
     * @param directory
     *          the directory
     * @throws IOException
     *          the IO exception
     */
    public static void writeStainsAsCsvData(
            String[] strains,
            File directory) throws IOException
    {
        File csvFile = new File(directory, STRAIN_NAMES_FILENAME);
        FlatFileWriter ffw = new FlatFileWriter(
                new FileWriter(csvFile),
                CommonFlatFileFormat.CSV_UNIX);
        ffw.writeRow(strains);
        ffw.flush();
        ffw.close();
    }
    
    /**
     * Bit-pack the input data and write a bunch of output files to the given
     * directory
     * @param parser
     *          the parser to use
     * @param inputCsvFile
     *          the input csv file
     * @param outputDirectory
     *          the output directory
     * @throws IOException
     *          if we run into trouble reading or writing data
     * @throws IllegalFormatException
     *          if an existing strain file in the directory has bad formatting
     */
    public static void writeHMMStatesAsBinaryData(
            HiddenMarkovModelStateParser parser,
            File inputCsvFile,
            File outputDirectory)
    throws IOException, IllegalFormatException
    {
        if(!outputDirectory.exists())
        {
            if(!outputDirectory.mkdirs())
            {
                throw new IOException(
                        "Failed to create directory: " +
                        outputDirectory.getAbsolutePath());
            }
        }
        else if(!outputDirectory.isDirectory())
        {
            throw new IOException(
                    outputDirectory.getAbsolutePath() +
                    " is not a valid directory");
        }
        
//        // read the states from CSV
//        BufferedReader hmmFileReader =
//            new BufferedReader(new FileReader(inputCsvFile));
//        Map<Integer, List<MultiPartitionedInterval>> states = parser.parseHMMStatesFromReader(
//                hmmFileReader,
//                null);
//        hmmFileReader.close();
//        
//        // read the strains from CSV
//        BufferedInputStream strainNameInStream =
//            new BufferedInputStream(new FileInputStream(inputCsvFile));
//        Set<String> strainSet = parser.parseAvailableStrainsFromStream(
//                strainNameInStream);
//        strainNameInStream.close();
//        String[] strains = strainSet.toArray(new String[strainSet.size()]);
//        Arrays.sort(strains);
//        
//        // write the states to binary
//        for(Entry<Integer, List<MultiPartitionedInterval>> hmmStatesEntry:
//            states.entrySet())
//        {
//            Integer currChr = hmmStatesEntry.getKey();
//            File chrDir = new File(
//                    outputDirectory,
//                    CHROMOSOME_DIR_PREFIX + currChr.toString());
//            chrDir.mkdir();
//            
//            List<MultiPartitionedInterval> currHMMStates =
//                hmmStatesEntry.getValue();
//            for(int i = 0; i < strains.length; i++)
//            {
//                String strain = strains[i];
//                try
//                {
//                    File strainFile = new File(
//                            chrDir,
//                            URL_CODEC.encode(strain) +
//                            HMM_STATE_STREAM_FILTER.getEndingString());
//                    
//                    if(strainFile.exists())
//                    {
//                        throw new IOException(
//                                "Refusing to overwrite file: " +
//                                strainFile.getAbsolutePath() +
//                                ". Please manually delete files or choose an empty " +
//                                "export directory.");
//                    }
//                    else
//                    {
//                        ObjectOutputStream oos = new ObjectOutputStream(
//                                new BufferedOutputStream(new FileOutputStream(
//                                        strainFile)));
//                        
//                        for(MultiPartitionedInterval currHMMInterval: currHMMStates)
//                        {
//                            oos.writeShort(currHMMInterval.get)
//                        }
//                        
//                        oos.flush();
//                        oos.close();
//                    }
//                }
//                catch(EncoderException ex)
//                {
//                    // hide it and throw it. we shouldn't run into this
//                    // exception but i want it to hurt if we do
//                    throw new RuntimeException(ex);
//                }
//            }
//        }
////////////////////////////////////////////////////////////////////////////////
        // read the strains from CSV
        BufferedInputStream strainNameInStream =
            new BufferedInputStream(new FileInputStream(inputCsvFile));
        Set<String> strainSet = parser.parseAvailableStrainsFromStream(
                strainNameInStream);
        strainNameInStream.close();
        String[] strains = strainSet.toArray(new String[strainSet.size()]);
        Arrays.sort(strains);
        
        File strainNamesOutputFile = new File(
                outputDirectory,
                STRAIN_NAMES_FILENAME);
        if(strainNamesOutputFile.exists())
        {
            // confirm that the strain names match up
            FlatFileReader ffr = new FlatFileReader(
                    new FileReader(strainNamesOutputFile),
                    CommonFlatFileFormat.CSV_UNIX);
            
            Set<String> existingStrains = new HashSet<String>(Arrays.asList(
                    ffr.readRow()));
            if(!strainSet.equals(existingStrains))
            {
                throw new IOException(
                        "Strain sets  should match across all chromosome " +
                        "files, but [" +
                        SequenceUtilities.toString(strainSet, ", ") +
                        "] does not match [" +
                        SequenceUtilities.toString(existingStrains, ", ") +
                        "].");
            }
        }
        else
        {
            // write the strains to CSV
            FlatFileWriter ffw = new FlatFileWriter(
                    new FileWriter(strainNamesOutputFile),
                    CommonFlatFileFormat.CSV_UNIX);
            ffw.writeRow(strains);
            ffw.flush();
            ffw.close();
        }
    }
    
//    /**
//     * Read in haplotype data using the given dir and chromosome number
//     * @param directory
//     *          the dir
//     * @param chromosomeNumber
//     *          the chromosome number
//     * @return
//     *          the list of intervals
//     * @throws FileNotFoundException
//     *          if the file doesn't exist (ie the chromosome isn't valid from
//     *          {@link #getAvailableBinaryDataChromosomes(File)}
//     * @throws IOException
//     *          if there is an IO failure
//     */
//    @SuppressWarnings("unchecked")
//    public static List<MultiPartitionedInterval> readBinaryMultiGroupHaplotypes(
//            File directory,
//            int chromosomeNumber) throws FileNotFoundException, IOException
//    {
//        File objectFile = new File(
//                directory,
//                CHROMOSOME_DIR_PREFIX + chromosomeNumber + CHROMOSOME_FILE_SUFFIX);
//        ObjectInputStream ois = new ObjectInputStream(new BufferedInputStream(
//                new FileInputStream(objectFile)));
//        
//        try
//        {
//            List<MultiPartitionedInterval> haplotypes =
//                (List<MultiPartitionedInterval>)ois.readObject();
//            
//            return haplotypes;
//        }
//        catch(ClassNotFoundException ex)
//        {
//            // this should never happen! throw it as a runtime
//            throw new RuntimeException(ex);
//        }
//        finally
//        {
//            ois.close();
//        }
//    }
    
    /**
     * Figure out which chromosome numbers are available in the given
     * directory
     * @param directory
     *          the directory
     * @return
     *          the chromosomes
     */
    public static int[] getAvailableBinaryDataChromosomes(File directory)
    {
        String[] fileNames = directory.list();
        List<Integer> chromosomes = new ArrayList<Integer>();
        for(String currFileName: fileNames)
        {
            if(currFileName.startsWith(CHROMOSOME_DIR_PREFIX))
            {
                String chrNumString = currFileName.substring(
                        CHROMOSOME_DIR_PREFIX.length());
                chromosomes.add(Integer.parseInt(chrNumString));
            }
        }
        
        int[] chrNums = SequenceUtilities.toIntArray(chromosomes);
        Arrays.sort(chrNums);
        return chrNums;
    }
}
