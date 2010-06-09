/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * Permission is hereby granted, free of charge, to any person obtaining  a copy
 * of this software and associated documentation files (the  "Software"), to
 * deal in the Software without restriction, including  without limitation the
 * rights to use, copy, modify, merge, publish,  distribute, sublicense, and/or
 * sell copies of the Software, and to  permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be  included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF  MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY  CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,  TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE  SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package org.jax.haplotype.io;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFileChooser;

import org.apache.commons.codec.EncoderException;
import org.apache.commons.codec.net.URLCodec;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.data.StreamingBinaryChromosomeDataSource;

/**
 * Some utility functions for SNP streams
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SnpStreamUtil
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            SnpStreamUtil.class.getName());
    
    private static final URLCodec URL_CODEC = new URLCodec();
    
    private static final String CHROMOSOME_DIR_PREFIX = "chr";
    
    /**
     * Write the strain chromosome to the snp stream
     * @param referenceStrainChromosome
     *          the reference
     * @param strainChromosomeToWrite
     *          the chromosome to write
     * @param outputStream
     *          the file to write to
     * @param streamDirection
     *          the stream direction to use
     * @throws IOException
     *          the IO stream
     */
    public static void writeBinarySnps(
            StrainChromosome referenceStrainChromosome,
            StrainChromosome strainChromosomeToWrite,
            OutputStream outputStream,
            StreamDirection streamDirection)
    throws IOException
    {
        final SnpInputStream snpInputStream;
        if(streamDirection == StreamDirection.FORWARD)
        {
            snpInputStream = new ForwardStrainChromosomeSnpInputStream(
                    referenceStrainChromosome,
                    strainChromosomeToWrite);
        }
        else
        {
            assert streamDirection == StreamDirection.REVERSE;
            snpInputStream = new ReverseStrainChromosomeSnpInputStream(
                    referenceStrainChromosome,
                    strainChromosomeToWrite);
        }
        
        SnpOutputStream snpOutputStream = new BinarySnpOutputStream(
                snpInputStream.getReadDirection(),
                outputStream,
                snpInputStream.getSnpCount());
        while(snpInputStream.hasNextSnp())
        {
            snpOutputStream.writeSnp(snpInputStream.getNextSnp());
        }
    }
    
    /**
     * Write the genotype data from the input csv file to a binary format
     * @param parser
     *          the parser to use
     * @param inputCsvFile
     *          the input
     * @param outputDirectory
     *          the output dir
     * @param strainsToWrite
     *          the strains
     * @param referenceStrainChromosome
     *          the reference strain (can be null)
     * @return
     *          the reference strain (one is selected if the given reference
     *          strain is null)
     * @throws IOException
     *          if we get one
     */
    private static StrainChromosome writeBinaryStrains(
            GenotypeParser parser,
            File inputCsvFile,
            File outputDirectory,
            Set<String> strainsToWrite,
            StrainChromosome referenceStrainChromosome)
    throws IOException
    {
        System.out.println(
                "Extracting " + strainsToWrite.size() +
                " strains from: " +
                inputCsvFile.getName());
        FileInputStream fis = new FileInputStream(inputCsvFile);
        Set<StrainChromosome> chromoSet = parser.parseGenotypeFromStream(
                fis,
                true,
                strainsToWrite);
        fis.close();
        
        Iterator<StrainChromosome> chromoIter = chromoSet.iterator();
        if(referenceStrainChromosome == null && chromoIter.hasNext())
        {
            referenceStrainChromosome = chromoIter.next();
            System.out.println(
                    "Reference Strain: " +
                    referenceStrainChromosome.getStrainName());
            
            // reset the iterator since we want to make sure we also write
            // out the reference strain
            chromoIter = chromoSet.iterator();
        }
        
        try
        {
            while(chromoIter.hasNext())
            {
                StrainChromosome nextChromo = chromoIter.next();
                File nextOutputDir = new File(
                        outputDirectory,
                        CHROMOSOME_DIR_PREFIX + nextChromo.getChromosomeNumber());
                nextOutputDir.mkdir();
                
                // Write the forward stream
                File nextOutputFile = new File(
                        nextOutputDir,
                        URL_CODEC.encode(nextChromo.getStrainName()) +
                        StreamingBinaryChromosomeDataSource.SNP_STREAM_FILTER.getEndingString());
                if(!nextOutputFile.exists())
                {
                    OutputStream nextOutputStream = new BufferedOutputStream(
                            new FileOutputStream(nextOutputFile));
                    SnpStreamUtil.writeBinarySnps(
                            referenceStrainChromosome,
                            nextChromo,
                            nextOutputStream,
                            StreamDirection.FORWARD);
                    
                    nextOutputStream.close();
                }
                else
                {
                    throw new IOException(
                            "refusing to overwrite " + nextOutputFile.getAbsolutePath());
                }
                
                // Write the reverse stream
                File nextReverseOutputFile = new File(
                        nextOutputDir,
                        URL_CODEC.encode(nextChromo.getStrainName()) +
                        StreamingBinaryChromosomeDataSource.REVERSE_SNP_STREAM_FILTER.getEndingString());
                if(!nextReverseOutputFile.exists())
                {
                    OutputStream nextReverseOutputStream = new BufferedOutputStream(
                            new FileOutputStream(nextReverseOutputFile));
                    SnpStreamUtil.writeBinarySnps(
                            referenceStrainChromosome,
                            nextChromo,
                            nextReverseOutputStream,
                            StreamDirection.REVERSE);
                    
                    nextReverseOutputStream.close();
                }
                else
                {
                    throw new IOException(
                            "refusing to overwrite " + nextOutputFile.getAbsolutePath());
                }
            }
        }
        catch(EncoderException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to encode strain name",
                    ex);
            throw new IOException(ex.getMessage());
        }
        
        return referenceStrainChromosome;
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
     */
    public static void writeBinaryChromosomeData(
            GenotypeParser parser,
            File inputCsvFile,
            File outputDirectory)
    throws IOException
    {
        FileInputStream fis = new FileInputStream(inputCsvFile);
        Set<String> strains = parser.parseAvailableStrainsFromStream(fis);
        fis.close();
        System.out.println("Strain count: " + strains.size());
        
        Set<String> strainSubset = new HashSet<String>();
        
        StrainChromosome referenceStrainChromosome = null;
        for(String strainName: strains)
        {
            strainSubset.add(strainName);
            if(strainSubset.size() == 20)
            {
                referenceStrainChromosome = writeBinaryStrains(
                        parser,
                        inputCsvFile,
                        outputDirectory,
                        strainSubset,
                        referenceStrainChromosome);
                strainSubset.clear();
            }
        }
        
        if(!strainSubset.isEmpty())
        {
            referenceStrainChromosome = writeBinaryStrains(
                    parser,
                    inputCsvFile,
                    outputDirectory,
                    strainSubset,
                    referenceStrainChromosome);
            strainSubset.clear();
        }
        
        writeBinarySnpPositions(
                referenceStrainChromosome,
                outputDirectory);
    }
    
    /**
     * Write the snp positions as a binary file using the given base directory
     * @param chromosome
     *          the chromosome to write
     * @param outputDirectory
     *          the base output directory
     * @throws IOException
     *          if the write operation fails
     */
    private static void writeBinarySnpPositions(
            StrainChromosome chromosome,
            File outputDirectory) throws IOException
    {
        File chrOutputDir = new File(
                outputDirectory,
                CHROMOSOME_DIR_PREFIX + chromosome.getChromosomeNumber());
        chrOutputDir.mkdir();
        
        File outputFile = new File(
                chrOutputDir,
                StreamingBinaryChromosomeDataSource.SNP_POSITION_FILE_NAME);
        if(outputFile.exists())
        {
            throw new IOException(
                    "refusing to overwrite existing file: " +
                    outputFile.getAbsolutePath());
        }
        DataOutputStream output = new DataOutputStream(
                new BufferedOutputStream(new FileOutputStream(outputFile)));
        writeBinarySnpPositions(
                chromosome,
                output,
                StreamDirection.FORWARD);
        output.flush();
        output.close();
        
        File reverseOutputFile = new File(
                chrOutputDir,
                StreamingBinaryChromosomeDataSource.REVERSE_SNP_POSITION_FILE_NAME);
        if(reverseOutputFile.exists())
        {
            throw new IOException(
                    "refusing to overwrite existing file: " +
                    reverseOutputFile.getAbsolutePath());
        }
        DataOutputStream reverseOutput = new DataOutputStream(
                new BufferedOutputStream(new FileOutputStream(reverseOutputFile)));
        writeBinarySnpPositions(
                chromosome,
                reverseOutput,
                StreamDirection.REVERSE);
        reverseOutput.flush();
        reverseOutput.close();
    }

    /**
     * Write the snp positions as a binary file using the given base directory
     * @param chromosome
     *          the chromosome to write
     * @throws IOException
     *          if the write operation fails
     */
    private static void writeBinarySnpPositions(
            StrainChromosome chromosome,
            DataOutputStream dataOutputStream,
            StreamDirection streamDirection)
    throws IOException
    {
        SingleNucleotidePolymorphism[] snps =
            chromosome.getSingleNucleotidePolymorphisms();
        dataOutputStream.write(StreamDirection.streamDirectionToByte(
                streamDirection));
        dataOutputStream.writeInt(chromosome.getChromosomeNumber());
        dataOutputStream.writeLong(snps[0].getPositionInBasePairs());
        dataOutputStream.writeLong(
                1L + snps[snps.length - 1].getPositionInBasePairs() -
                snps[0].getPositionInBasePairs());
        dataOutputStream.writeLong(
                snps.length);
        if(streamDirection == StreamDirection.FORWARD)
        {
            for(int i = 0; i < snps.length; i++)
            {
                dataOutputStream.writeLong(snps[i].getPositionInBasePairs());
            }
        }
        else
        {
            assert streamDirection == StreamDirection.REVERSE;
            for(int i = snps.length - 1; i >= 0; i--)
            {
                dataOutputStream.writeLong(snps[i].getPositionInBasePairs());
            }
        }
    }

    /**
     * Create chromosome data sources from the given binary data directory
     * @param dataDirectory
     *          the data dir that contains chromosome data sources
     * @return
     *          the mapping from chromosome number to chromosome data source
     */
    public static Map<Integer, StreamingBinaryChromosomeDataSource> getBinaryChromosomeDataSources(
            File dataDirectory)
    {
        File[] allFiles = dataDirectory.listFiles();
        
        if(allFiles == null)
        {
            throw new NullPointerException(
                    "Cannot list chromosome directories in: " +
                    dataDirectory.getAbsolutePath());
        }
        
        Map<Integer, StreamingBinaryChromosomeDataSource> chromosomeDataSources =
            new HashMap<Integer, StreamingBinaryChromosomeDataSource>();
        for(File file: allFiles)
        {
            if(file.isDirectory() && file.getName().startsWith(CHROMOSOME_DIR_PREFIX))
            {
                String chromosomeNumberString =
                    file.getName().substring(CHROMOSOME_DIR_PREFIX.length());
                try
                {
                    Integer chromosomeNumber = Integer.valueOf(
                            chromosomeNumberString);
                    StreamingBinaryChromosomeDataSource chromosomeDataSource =
                        new StreamingBinaryChromosomeDataSource(
                                file,
                                chromosomeNumber.intValue());
                    chromosomeDataSources.put(
                            chromosomeNumber,
                            chromosomeDataSource);
                }
                catch(NumberFormatException ex)
                {
                    LOG.log(Level.WARNING,
                            "Directory doesn't match chromosome format: " +
                            file.getAbsolutePath() + ". Ignoring...",
                            ex);
                }
            }
        }
        return chromosomeDataSources;
    }

    /**
     * A main for snp data conversion
     * @param args
     *          dont care
     * @throws IOException
     *          if we get one
     */
    public static void main(String[] args) throws IOException
    {
        JFileChooser inputFileChooser = new JFileChooser();
        inputFileChooser.setDialogTitle("Select CSV Chromosome Input Files");
        inputFileChooser.setMultiSelectionEnabled(true);
        inputFileChooser.setFileSelectionMode(
                JFileChooser.FILES_ONLY);
        int userSelection = inputFileChooser.showOpenDialog(null);
        if(userSelection == JFileChooser.APPROVE_OPTION)
        {
            JFileChooser outputDirectoryChooser = new JFileChooser();
            outputDirectoryChooser.setDialogTitle(
                    "Select an Output Directory");
            outputDirectoryChooser.setFileSelectionMode(
                    JFileChooser.DIRECTORIES_ONLY);
            int outputUserSelection =
                outputDirectoryChooser.showOpenDialog(null);
            
            if(outputUserSelection == JFileChooser.APPROVE_OPTION)
            {
                File selectedOutputDirectory = outputDirectoryChooser.getSelectedFile();
                File[] selectedInputFiles = inputFileChooser.getSelectedFiles();
                for(File selectedInputFile: selectedInputFiles)
                {
                    writeBinaryChromosomeData(
                            new GenotypeParser(),
                            selectedInputFile,
                            selectedOutputDirectory);
                }
            }
        }
        else
        {
            System.out.println("user doesn't want to open the file");
        }
    }
}
