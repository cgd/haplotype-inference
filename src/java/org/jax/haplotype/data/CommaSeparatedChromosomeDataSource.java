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

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.ForwardStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.ForwardStrainChromosomeSnpPositionInputStream;
import org.jax.haplotype.io.GenotypeParser;
import org.jax.haplotype.io.ReverseStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.ReverseStrainChromosomeSnpPositionInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SimpleSdpInputStream;
import org.jax.haplotype.io.SnpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.haplotype.io.StreamDirection;

/**
 * A data source for comma-separated chromosome data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommaSeparatedChromosomeDataSource implements ChromosomeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -3873922254783448637L;

    private static final Logger LOG = Logger.getLogger(
            CommaSeparatedChromosomeDataSource.class.getName());
    
    private final URL csvGenotypeUrl;
    
    private volatile int chromosomeNumber;

    private Set<String> cachedStrainNames;
    
    private Map<String, StrainChromosome> cachedGenotypeData;

    private final boolean cacheData;
    
    /**
     * Constructor with no data caching
     * @param csvGenotypeFile
     *          the comma seperated file that we get the data from
     * @param chromosomeNumber
     *          the chromosome number (-1 indicates that the chromosome
     *          should be read from the data)
     */
    public CommaSeparatedChromosomeDataSource(
            File csvGenotypeFile,
            int chromosomeNumber)
    {
        this(csvGenotypeFile, chromosomeNumber, false);
    }
    
    /**
     * Constructor
     * @param csvGenotypeFile
     *          the comma seperated file that we get the data from
     * @param chromosomeNumber
     *          the chromosome number (-1 indicates that the chromosome
     *          should be read from the data)
     * @param cacheData
     *          if true cache all of the data that is uploaded
     */
    public CommaSeparatedChromosomeDataSource(
            File csvGenotypeFile,
            int chromosomeNumber,
            boolean cacheData)
    {
        try
        {
            this.csvGenotypeUrl = csvGenotypeFile.toURI().toURL();
            this.chromosomeNumber = chromosomeNumber;
            this.cacheData = cacheData;
            this.cachedGenotypeData = null;
            this.cachedStrainNames = null;
        }
        catch(MalformedURLException ex)
        {
            throw new RuntimeException(ex);
        }
    }
    
    /**
     * Constructor with no data caching
     * @param csvGenotypeUrl
     *          the comma seperated file that we get the data from
     * @param chromosomeNumber
     *          the chromosome number (-1 indicates that the chromosome
     *          should be read from the data)
     */
    public CommaSeparatedChromosomeDataSource(
            URL csvGenotypeUrl,
            int chromosomeNumber)
    {
        this(csvGenotypeUrl, chromosomeNumber, false);
    }
    
    /**
     * Constructor
     * @param csvGenotypeUrl
     *          the comma seperated file that we get the data from
     * @param chromosomeNumber
     *          the chromosome number (-1 indicates that the chromosome
     *          should be read from the data)
     * @param cacheData
     *          if true cache all of the data that is uploaded
     */
    public CommaSeparatedChromosomeDataSource(
            URL csvGenotypeUrl,
            int chromosomeNumber,
            boolean cacheData)
    {
        this.csvGenotypeUrl = csvGenotypeUrl;
        this.chromosomeNumber = chromosomeNumber;
        this.cacheData = cacheData;
        this.cachedGenotypeData = null;
        this.cachedStrainNames = null;
    }
    
    /**
     * Getter for the backing CSV url
     * @return the csvGenotypeUrl
     */
    public URL getCsvGenotypeUrl()
    {
        return this.csvGenotypeUrl;
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getAvailableStrains()
    {
        try
        {
            if(this.cacheData)
            {
                if(this.cachedStrainNames != null)
                {
                    return new HashSet<String>(this.cachedStrainNames);
                }
                else
                {
                    this.cachedStrainNames = this.getAvailableStrains(this.csvGenotypeUrl);
                    return new HashSet<String>(this.cachedStrainNames);
                }
            }
            else
            {
                return this.getAvailableStrains(this.csvGenotypeUrl);
            }
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to get available strains",
                    ex);
            return null;
        }
    }

    /**
     * {@inheritDoc}
     */
    public Set<StrainChromosome> getGenotypeData(Set<String> strainsToParse)
    {
        try
        {
            final Set<String> newStrainsToParse;
            
            // remove any cached values that we already know about
            if(this.cacheData && this.cachedGenotypeData != null)
            {
                newStrainsToParse = new HashSet<String>(strainsToParse);
                newStrainsToParse.removeAll(this.cachedGenotypeData.keySet());
            }
            else
            {
                newStrainsToParse = strainsToParse;
            }
            
            Set<StrainChromosome> newGenotypeData;
            if(newStrainsToParse.isEmpty())
            {
                newGenotypeData = Collections.emptySet();
            }
            else
            {
                newGenotypeData =
                    this.getGenotypeData(this.csvGenotypeUrl, newStrainsToParse);
            }
            
            if(this.cacheData)
            {
                if(this.cachedGenotypeData == null)
                {
                    this.cachedGenotypeData = Collections.synchronizedMap(
                            new HashMap<String, StrainChromosome>());
                }
                
                for(StrainChromosome newChromosome: newGenotypeData)
                {
                    this.cachedGenotypeData.put(
                            newChromosome.getStrainName(),
                            newChromosome);
                }
                
                Set<StrainChromosome> genotypeData = new HashSet<StrainChromosome>(
                        strainsToParse.size());
                for(String strainName: strainsToParse)
                {
                    genotypeData.add(this.cachedGenotypeData.get(strainName));
                }
                return genotypeData;
            }
            else
            {
                return newGenotypeData;
            }
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to get genotype data",
                    ex);
            return null;
        }
    }
    
    /**
     * Getter for the genotype data that we're going to infer the snps from
     * @param genotypeDataSource
     *          the genotype data source
     * @param strainsToParse
     *          the strains to parse
     * @return
     *          the data
     */
    private Set<StrainChromosome> getGenotypeData(
            URL genotypeDataSource,
            Set<String> strainsToParse)
    throws IOException
    {
        InputStream is = genotypeDataSource.openStream();
        GenotypeParser parser =
            new GenotypeParser();
        return parser.parseGenotypeFromStream(
                is,
                true,
                strainsToParse);
    }
    
    /**
     * Extract the available strains from a genotype file
     * @param genotypeDataSource
     *          the URL
     * @return
     *          the strain name set
     */
    private Set<String> getAvailableStrains(
            URL genotypeDataSource)
    throws IOException
    {
        InputStream is = genotypeDataSource.openStream();
        GenotypeParser parser =
            new GenotypeParser();
        return parser.parseAvailableStrainsFromStream(is);
    }
    
    /**
     * Get the available strains in the order that they appear in the CSV
     * file
     * @return
     *          the strains
     */
    public List<String> getAvailableStrainsInOrder()
    {
        try
        {
            return this.getAvailableStrainsInOrder(this.csvGenotypeUrl);
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to get available strains",
                    ex);
            return null;
        }
    }
    
    /**
     * Get indices for the strain set
     * @param strains
     *          the strains
     * @return
     *          the indices
     */
    public List<Integer> getStrainIndices(
            Collection<String> strains)
    {
        return CommaSeparatedChromosomeDataSource.getStrainIndices(
                strains,
                this.getAvailableStrainsInOrder());
    }

    /**
     * Helper function to
     * get a strain index list for the strain set using the given ordering
     * @param strains
     *          the strains to find indices for
     * @param orderedStrains
     *          strains in index order
     * @return
     *          the indices
     */
    public static List<Integer> getStrainIndices(
            Collection<String> strains,
            List<String> orderedStrains)
    {
        List<Integer> indices = new ArrayList<Integer>();
        for(String strain: strains)
        {
            int index = orderedStrains.indexOf(strain);
            if(index != -1)
            {
                indices.add(index);
            }
        }
        Collections.sort(indices);
        return indices;
    }
    
    /**
     * Get the available strains in the order that they appear in the CSV
     * file
     * @param genotypeDataSource
     *          the file to read the strain names from
     * @return
     *          the strains
     */
    private List<String> getAvailableStrainsInOrder(URL genotypeDataSource)
    throws IOException
    {
        InputStream is = genotypeDataSource.openStream();
        GenotypeParser parser =
            new GenotypeParser();
        return parser.parseAvailableStrainsInOrderFromStream(is);
    }
    
    /**
     * {@inheritDoc}
     */
    public int getChromosomeNumber()
    {
        if(this.chromosomeNumber == -1)
        {
            this.chromosomeNumber = this.getAnyChromosome().getChromosomeNumber();
            
            if(this.chromosomeNumber == -1)
            {
                throw new IllegalStateException(
                        "Failed to read chromosome number from URL");
            }
        }
        
        return this.chromosomeNumber;
    }
    
    /**
     * {@inheritDoc}
     */
    public long getDataExtentInBasePairs()
    {
        // TODO we really need caching here
        Set<StrainChromosome> genoData = this.getGenotypeData(
                Collections.singleton(this.getAvailableStrains().iterator().next()));
        StrainChromosome firstChromosome = genoData.iterator().next();
        SingleNucleotidePolymorphism[] snps =
            firstChromosome.getSingleNucleotidePolymorphisms();
        return snps[snps.length - 1].getPositionInBasePairs() -
               snps[0].getPositionInBasePairs();
    }

    /**
     * {@inheritDoc}
     */
    public long getDataStartInBasePairs()
    {
        // TODO we really need caching here
        Set<StrainChromosome> genoData = this.getGenotypeData(
                Collections.singleton(this.getAvailableStrains().iterator().next()));
        StrainChromosome firstChromosome = genoData.iterator().next();
        return firstChromosome.getSingleNucleotidePolymorphisms()[0].getPositionInBasePairs();
    }
    
    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String referenceStrainName,
            String[] comparisonStrainNames)
    {
        Set<String> strainsToParse = new HashSet<String>();
        for(int i = 0; i < comparisonStrainNames.length; i++)
        {
            strainsToParse.add(comparisonStrainNames[i]);
        }
        strainsToParse.add(referenceStrainName);
        
        SnpInputStream[] snpInputStreams =
            new SnpInputStream[comparisonStrainNames.length];
        Set<StrainChromosome> genoData = this.getGenotypeData(strainsToParse);
        Map<String, StrainChromosome> genoMap =
            new HashMap<String, StrainChromosome>(genoData.size());
        for(StrainChromosome strainChromosome: genoData)
        {
            genoMap.put(strainChromosome.getStrainName(), strainChromosome);
        }
        
        switch(streamDirection)
        {
            case FORWARD:
            {
                for(int i = 0; i < snpInputStreams.length; i++)
                {
                    snpInputStreams[i] = new ForwardStrainChromosomeSnpInputStream(
                            genoMap.get(referenceStrainName),
                            genoMap.get(comparisonStrainNames[i]));
                }
            }
            break;
            
            case REVERSE:
            {
                for(int i = 0; i < snpInputStreams.length; i++)
                {
                    snpInputStreams[i] = new ReverseStrainChromosomeSnpInputStream(
                            genoMap.get(referenceStrainName),
                            genoMap.get(comparisonStrainNames[i]));
                }
            }
            break;
            
            default:
            {
                throw new IllegalArgumentException(
                        "unknown stream direction: " + streamDirection);
            }
        }
        
        try
        {
            return new SimpleSdpInputStream(
                    comparisonStrainNames,
                    snpInputStreams);
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE,
                    "Failed to construct a simple SDP stream",
                    ex);
            return null;
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
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String[] strainNames)
    {
        return this.getSdpInputStream(
                streamDirection,
                strainNames[0],
                strainNames);
    }
    
    /**
     * {@inheritDoc}
     */
    public SdpInputStream getSdpInputStream(String[] strainNames)
    {
        return this.getSdpInputStream(strainNames[0], strainNames);
    }
    
    /**
     * {@inheritDoc}
     */
    public SnpPositionInputStream getSnpPositionInputStream(
            StreamDirection streamDirection)
    {
        StrainChromosome anyChromosome = this.getAnyChromosome();
        
        switch(streamDirection)
        {
            case FORWARD:
            {
                return new ForwardStrainChromosomeSnpPositionInputStream(
                        anyChromosome);
            }
            
            case REVERSE:
            {
                return new ReverseStrainChromosomeSnpPositionInputStream(
                        anyChromosome);
            }
            
            default:
            {
                throw new IllegalArgumentException(
                        "unknown stream direction: " + streamDirection);
            }
        }
    }
    
    /**
     * Get any of the chromosome data
     * @return
     *          any chromosome
     */
    private StrainChromosome getAnyChromosome()
    {
        if(this.cachedGenotypeData == null || this.cachedGenotypeData.isEmpty())
        {
            List<String> orderedStrains = this.getAvailableStrainsInOrder();
            Set<StrainChromosome> anyGenoData = this.getGenotypeData(
                    Collections.singleton(orderedStrains.get(0)));
            return anyGenoData.iterator().next();
        }
        else
        {
            return this.cachedGenotypeData.values().iterator().next();
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
