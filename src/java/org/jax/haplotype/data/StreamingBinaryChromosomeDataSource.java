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
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import org.apache.commons.codec.DecoderException;
import org.apache.commons.codec.net.URLCodec;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.BinarySnpInputStream;
import org.jax.haplotype.io.BinarySnpPositionInputStream;
import org.jax.haplotype.io.ReferenceNormalizedSdpInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SimpleSdpInputStream;
import org.jax.haplotype.io.SnpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.haplotype.io.StreamDirection;
import org.jax.util.io.FileExtensionFilter;

/**
 * A {@link ChromosomeDataSource} that works off of streaming binary data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class StreamingBinaryChromosomeDataSource implements ChromosomeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -215758818785132881L;

    /**
     * File filter used for the forward SNP stream
     */
    public static final FileExtensionFilter SNP_STREAM_FILTER =
        new FileExtensionFilter("snp");
    
    /**
     * File filter used for the forward SNP stream
     */
    public static final FileExtensionFilter REVERSE_SNP_STREAM_FILTER =
        new FileExtensionFilter("rsnp");
    
    /**
     * File used for the forward SNP position stream
     */
    public static final String SNP_POSITION_FILE_NAME =
        "snp.pos";
    
    /**
     * File used for the reverse SNP position stream
     */
    public static final String REVERSE_SNP_POSITION_FILE_NAME =
        "snp.rpos";
    
    private static final URLCodec URL_CODEC = new URLCodec();

    private final int chromosomeNumber;

    private final File dataDirectory;
    
    private volatile Set<String> persistentStrainsToAcceptFilter = null;
    
    /**
     * Constructor
     * @param dataDirectory
     *          the data directory for this chromosome data source
     * @param chromosomeNumber
     *          the chromosome number for this chromosome data source
     */
    public StreamingBinaryChromosomeDataSource(
            File dataDirectory,
            int chromosomeNumber)
    {
        this.chromosomeNumber = chromosomeNumber;
        this.dataDirectory = dataDirectory;
    }
    
    /**
     * Getter for the persistent strain filter
     * @return the persistent filter
     */
    public Set<String> getPersistentStrainsToAcceptFilter()
    {
        return this.persistentStrainsToAcceptFilter;
    }
    
    /**
     * Setter for the persistent filter
     * @param persistentStrainsToAcceptFilter
     *          the persistent filter
     */
    public void setPersistentStrainsToAcceptFilter(
            Set<String> persistentStrainsToAcceptFilter)
    {
        this.persistentStrainsToAcceptFilter = persistentStrainsToAcceptFilter;
    }
    
    /**
     * This is the same as {@link #getAvailableStrains()} except that the
     * {@link #getPersistentStrainsToAcceptFilter()} is not applied to the
     * strain names
     * @return
     *          all available strains whether they're in the filter list
     *          or not
     */
    public Set<String> getUnfilteredStrains()
    {
        Map<String, File> strainNameToFileMap = this.getStrainNameToSnpFileMap(
                SNP_STREAM_FILTER,
                null);
        
        return strainNameToFileMap.keySet();
    }
    
    /**
     * {@inheritDoc}
     */
    public Set<String> getAvailableStrains()
    {
        Map<String, File> strainNameToFileMap = this.getStrainNameToSnpFileMap(
                SNP_STREAM_FILTER,
                this.persistentStrainsToAcceptFilter);
        
        return strainNameToFileMap.keySet();
    }
    
    private Map<String, File> getStrainNameToSnpFileMap(
            FileExtensionFilter fileFilter,
            Set<String> strainNameFilter)
    {
        try
        {
            File[] strainFiles = this.dataDirectory.listFiles(fileFilter);
            Map<String, File> strainFileMap = new HashMap<String, File>(
                    strainFiles.length);
            for(File file: strainFiles)
            {
                String name = file.getName();
                String nameSansExtension = name.substring(
                        0,
                        name.length() - fileFilter.getEndingString().length());
                String strainName = StreamingBinaryChromosomeDataSource.URL_CODEC.decode(
                        nameSansExtension);
                
                if(strainNameFilter == null || strainNameFilter.contains(strainName))
                {
                    strainFileMap.put(strainName, file);
                }
            }
            
            return strainFileMap;
        }
        catch(DecoderException ex)
        {
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public int getChromosomeNumber()
    {
        return this.chromosomeNumber;
    }

    /**
     * {@inheritDoc}
     */
    public long getDataExtentInBasePairs()
    {
        return this.getSnpPositionInputStream().getExtentInBasePairs();
    }

    /**
     * {@inheritDoc}
     */
    public long getDataStartInBasePairs()
    {
        return this.getSnpPositionInputStream().getStartInBasePairs();
    }

    /**
     * {@inheritDoc}
     */
    public Set<StrainChromosome> getGenotypeData(Set<String> strainsToParse)
    {
        // TODO we need the observed call stream to be able to implement this
        throw new UnsupportedOperationException(
                "Not yet implemented");
    }
    
    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String[] strainNames)
    {
        try
        {
            Map<String, File> strainToFileMap =
                this.getStrainNameToSnpFileMap(streamDirection);
            
            SnpInputStream[] snpInputStreams = new SnpInputStream[strainNames.length];
            for(int i = 0; i < snpInputStreams.length; i++)
            {
                snpInputStreams[i] = new BinarySnpInputStream(
                        new BufferedInputStream(new FileInputStream(
                                strainToFileMap.get(strainNames[i]))));
            }
            
            return new SimpleSdpInputStream(
                    strainNames,
                    snpInputStreams);
        }
        catch(IOException ex)
        {
            throw new RuntimeException(ex);
        }
    }

    /**
     * Get the mapping from strain name to SNP file given the stream
     * direction that we're looking for
     * @param streamDirection
     *          the SNP stream direction
     * @return
     *          the strain name to SNP file mapping
     */
    private Map<String, File> getStrainNameToSnpFileMap(
            StreamDirection streamDirection)
    {
        switch(streamDirection)
        {
            case FORWARD:
            {
                return this.getStrainNameToSnpFileMap(
                        SNP_STREAM_FILTER,
                        this.persistentStrainsToAcceptFilter);
            }
            
            case REVERSE:
            {
                return this.getStrainNameToSnpFileMap(
                        REVERSE_SNP_STREAM_FILTER,
                        this.persistentStrainsToAcceptFilter);
            }
            
            default:
            {
                throw new IllegalArgumentException(
                        "Unknown stream direction: " +
                        streamDirection);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(String[] strainNames)
    {
        return this.getSdpInputStream(StreamDirection.FORWARD, strainNames);
    }
    
    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String referenceStrainName,
            String[] comparisonStrainNames)
    {
        try
        {
            Map<String, File> strainToFileMap =
                this.getStrainNameToSnpFileMap(streamDirection);
            
            SnpInputStream[] snpInputStreams = new SnpInputStream[comparisonStrainNames.length];
            for(int i = 0; i < snpInputStreams.length; i++)
            {
                snpInputStreams[i] = new BinarySnpInputStream(
                        new BufferedInputStream(new FileInputStream(
                                strainToFileMap.get(comparisonStrainNames[i]))));
            }
            
            SdpInputStream comparisonSdpStream = new SimpleSdpInputStream(
                    comparisonStrainNames,
                    snpInputStreams);
            SnpInputStream referenceSnpStream = new BinarySnpInputStream(
                    new BufferedInputStream(new FileInputStream(
                            strainToFileMap.get(referenceStrainName))));
            
            return new ReferenceNormalizedSdpInputStream(
                    referenceSnpStream,
                    comparisonSdpStream);
        }
        catch(IOException ex)
        {
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(
            String referenceStrainName,
            String[] comparisonStrainNames)
    {
        return this.getSdpInputStream(
                StreamDirection.FORWARD,
                referenceStrainName,
                comparisonStrainNames);
    }
    
    /**
     * {@inheritDoc}
     */
    public SnpPositionInputStream getSnpPositionInputStream(
            StreamDirection streamDirection)
    {
        File snpPositionsFile;
        switch(streamDirection)
        {
            case FORWARD:
            {
                snpPositionsFile = new File(
                        this.dataDirectory,
                        SNP_POSITION_FILE_NAME);
            }
            break;
            
            case REVERSE:
            {
                snpPositionsFile = new File(
                        this.dataDirectory,
                        REVERSE_SNP_POSITION_FILE_NAME);
            }
            break;
            
            default:
            {
                throw new IllegalArgumentException(
                        "Unknown stream direction: " +
                        streamDirection);
            }
        }
        
        try
        {
            return new BinarySnpPositionInputStream(new DataInputStream(
                    new BufferedInputStream(
                            new FileInputStream(snpPositionsFile))));
        }
        catch(IOException ex)
        {
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public SnpPositionInputStream getSnpPositionInputStream()
    {
        return this.getSnpPositionInputStream(StreamDirection.FORWARD);
    }
}
