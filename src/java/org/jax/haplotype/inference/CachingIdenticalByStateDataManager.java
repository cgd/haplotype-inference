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

package org.jax.haplotype.inference;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Pattern;

import org.apache.commons.codec.EncoderException;
import org.apache.commons.codec.net.URLCodec;
import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SnpIntervalList;
import org.jax.geneticutil.data.SnpIntervalListGroup;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.expressions.IntersectionIntervalExpression;
import org.jax.haplotype.expressions.IntervalCompliment;
import org.jax.haplotype.expressions.UnionIntervalExpression;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.util.io.IllegalFormatException;

/**
 * A caching implementation of the IBS finder
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CachingIdenticalByStateDataManager
{
    private static final Logger LOG = Logger.getLogger(
            CachingIdenticalByStateDataManager.class.getName());
    
    private static final CachingIdenticalByStateDataManager instance =
        new CachingIdenticalByStateDataManager(
                CachingGenomeDataManager.getInstance());
    
    private final CachingGenomeDataManager genomeDataManager;
    
    /**
     * used for encoding arbitrary strings as safe file names
     */
    protected static final URLCodec URL_CODEC = new URLCodec();
    
    /**
     * Used to join strain names
     */
    protected static final String CONCATINATION_STRING = "&";
    
    private static final String OR_REGEX = Pattern.quote(" or ");
    private static final String AND_REGEX = Pattern.quote(" and ");
    private static final String EQUALS_REGEX = Pattern.quote("=");
    private static final String NOT_EQUALS_REGEX = Pattern.quote("\u2260");
    
    /**
     * @param genomeDataManager
     */
    private CachingIdenticalByStateDataManager(
            CachingGenomeDataManager genomeDataManager)
    {
        this.genomeDataManager = genomeDataManager;
    }

    /**
     * Get the singleton instance
     * @return
     *          the instance
     */
    public static CachingIdenticalByStateDataManager getInstance()
    {
        return CachingIdenticalByStateDataManager.instance;
    }

    /**
     * Getter for the genome data manager used by this IBS manager
     * @return
     *          the genome data manager
     */
    public CachingGenomeDataManager getGenomeDataManager()
    {
        return this.genomeDataManager;
    }
    
    /**
     * Evaluate the strain comparison expression
     * @param genomeName
     *          the genome
     * @param strainComparisonExpression
     *          the expression
     * @param chromosomeNumber
     *          the chromosome number
     * @param minimumExtentInSnps
     *          the minimum extent in snps
     * @param minimumExtentInBasePairs
     *          the minimum extent in base pairs
     * @return
     *          the list of intervals where the given expression holds
     *          true for the given constraints
     */
    public SnpIntervalList evaluateStrainComparison(
            String genomeName,
            String strainComparisonExpression,
            int chromosomeNumber,
            int minimumExtentInSnps,
            long minimumExtentInBasePairs)
    {
        if(LOG.isLoggable(Level.FINE))
        {
            LOG.fine("evaluating: " + strainComparisonExpression);
        }
        
        try
        {
            // split into OR groups
            String[] orGroups = strainComparisonExpression.split(OR_REGEX);
            
            long intervalListStartBp = -1;
            long intervalListExtentBp = -1;
            
            BasePairInterval[] cumulativeOrIntervals = null;
            for(String orGroup: orGroups)
            {
                if(LOG.isLoggable(Level.FINE))
                {
                    LOG.fine("OR group: " + orGroup);
                }
                
                BasePairInterval[] cumulativeAndIntervals = null;
                
                // split into strain comparisons
                String[] strainComparisons = orGroup.split(AND_REGEX);
                for(String strainComparison: strainComparisons)
                {
                    if(LOG.isLoggable(Level.FINE))
                    {
                        LOG.fine("strain comparison: " + strainComparison);
                    }
                    
                    boolean calcNotEquals = false;
                    String[] strains = strainComparison.split(EQUALS_REGEX);
                    if(strains.length == 1)
                    {
                        calcNotEquals = true;
                        strains = strainComparison.split(NOT_EQUALS_REGEX);
                    }
                    
                    if(strains.length != 2)
                    {
                        throw new IllegalFormatException(
                                "bad strain comparison: " + strainComparison);
                    }
                    
                    String strain1 = strains[0];
                    String strain2 = strains[1];
                    if(LOG.isLoggable(Level.FINE))
                    {
                        if(calcNotEquals)
                        {
                            LOG.fine("calculating not equals for strains:");
                        }
                        else
                        {
                            LOG.fine("calculating equals for strains:");
                        }
                        
                        LOG.fine("strain 1: " + strain1);
                        LOG.fine("strain 2: " + strain2);
                    }
                    
                    SnpIntervalListGroup currGroup = this.findIdenticalByStateRegions(
                            genomeName,
                            strain1,
                            new String[] {strain2},
                            chromosomeNumber,
                            minimumExtentInSnps,
                            minimumExtentInBasePairs);
                    List<BasePairInterval> snpBlockList =
                        currGroup.getSnpBlocksMap().values().iterator().next();
                    intervalListStartBp = currGroup.getStartInBasePairs();
                    intervalListExtentBp = currGroup.getExtentInBasePairs();
                    
                    BasePairInterval[] currIntervals = snpBlockList.toArray(
                            new BasePairInterval[snpBlockList.size()]);
                    if(calcNotEquals)
                    {
                        currIntervals = IntervalCompliment.calculateCompliment(
                                currIntervals,
                                chromosomeNumber,
                                intervalListStartBp,
                                intervalListExtentBp);
                    }
                    
                    if(cumulativeAndIntervals == null)
                    {
                        cumulativeAndIntervals = currIntervals;
                    }
                    else
                    {
                        cumulativeAndIntervals = IntersectionIntervalExpression.intersect(
                                cumulativeAndIntervals,
                                currIntervals);
                    }
                }
                
                if(cumulativeOrIntervals == null)
                {
                    cumulativeOrIntervals = cumulativeAndIntervals;
                }
                else
                {
                    cumulativeOrIntervals = UnionIntervalExpression.calculateUnion(
                            cumulativeOrIntervals,
                            cumulativeAndIntervals);
                }
            }
            
            if(cumulativeOrIntervals == null)
            {
                LOG.severe(
                        "no intervals to calculate: " +
                        strainComparisonExpression);
                return null;
            }
            else
            {
                return new SnpIntervalList(
                        Arrays.asList(cumulativeOrIntervals),
                        intervalListStartBp,
                        intervalListExtentBp);
            }
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to evaluate strain comparison: " +
                    strainComparisonExpression,
                    ex);
            return null;
        }
    }
    
    /**
     * Find the IBS regions
     * @param genomeName
     *          the name of the genome that we want to use
     * @param referenceStrain
     *          the reference strain
     * @param comparisonStrains
     *          the comparison strains
     * @param chromosomeNumber
     *          the chromosome number (starting with 1)
     * @param minimumExtentInSnps
     *          the minimum IBD extent in SNPs
     * @param minimumExtentInBasePairs
     *          the minimum IBD extent in base pairs
     * @return
     *          the SNP blocks
     */
    public SnpIntervalListGroup findIdenticalByStateRegions(
            String genomeName,
            String referenceStrain,
            String[] comparisonStrains,
            int chromosomeNumber,
            int minimumExtentInSnps,
            long minimumExtentInBasePairs)
    {
        try
        {
            for(int i = 0; i < comparisonStrains.length; i++)
            {
                if(comparisonStrains[i] == null)
                {
                    throw new NullPointerException();
                }
            }
            
            ScanningIdenticalByStateFinder scanningIdenticalByStateFinder =
                new ScanningIdenticalByStateFinder();
            
            LOG.fine(
                    "calling caching findIdenticalByStateRegions. " +
                    "# comparison strains: " + comparisonStrains.length);
            
            Map<String, SnpIntervalList> snpBlockListsMap =
                new HashMap<String, SnpIntervalList>();
            synchronized(CachingIdenticalByStateDataManager.class)
            {
                Map<String, File> comparisonStrainsFileMap =
                    new HashMap<String, File>();
                for(int i = 0; i < comparisonStrains.length; i++)
                {
                    File uniqueFileCache = getCacheFile(
                            genomeName,
                            referenceStrain,
                            comparisonStrains[i],
                            chromosomeNumber,
                            minimumExtentInSnps,
                            minimumExtentInBasePairs);
                    
                    // if the cache doesn't exist, create it
                    if(uniqueFileCache.createNewFile())
                    {
                        uniqueFileCache.deleteOnExit();
                        comparisonStrainsFileMap.put(
                                comparisonStrains[i],
                                uniqueFileCache);
                    }
                    else
                    {
                        LOG.fine(
                                "loading cache: " +
                                uniqueFileCache.getAbsolutePath());
                        ObjectInputStream objectInput = new ObjectInputStream(
                                new BufferedInputStream(new FileInputStream(uniqueFileCache)));
                        snpBlockListsMap.put(
                                comparisonStrains[i],
                                (SnpIntervalList)objectInput.readObject());
                    }
                }
                
                if(!comparisonStrainsFileMap.isEmpty())
                {
                    GenomeDataSource genomeDataSource = this.genomeDataManager.getGenomeDataMap().get(
                            genomeName);
                    ChromosomeDataSource chromosomeDataSource =
                        genomeDataSource.getChromosomeDataSources().get(
                                chromosomeNumber);
                    
                    String[] comparisonStrainsToCalculate =
                        comparisonStrainsFileMap.keySet().toArray(
                                new String[comparisonStrainsFileMap.size()]);
                    for(int i = 0; i < comparisonStrainsToCalculate.length; i++)
                    {
                        if(comparisonStrainsToCalculate[i] == null)
                        {
                            throw new NullPointerException();
                        }
                    }
                    SdpInputStream sdpStream = chromosomeDataSource.getSdpInputStream(
                            referenceStrain,
                            comparisonStrainsToCalculate);
                    SnpPositionInputStream snpPositionStream =
                        chromosomeDataSource.getSnpPositionInputStream();
                    SnpIntervalListGroup newIbsRegions = scanningIdenticalByStateFinder.findIdenticalByStateRegions(
                            sdpStream,
                            snpPositionStream,
                            minimumExtentInSnps,
                            minimumExtentInBasePairs);
                    Map<String, List<BasePairInterval>> newIbsIntervalLists =
                        newIbsRegions.getSnpBlocksMap();
                    for(Entry<String, List<BasePairInterval>> intervalListEntry:
                        newIbsIntervalLists.entrySet())
                    {
                        System.out.println("Strain Key: " + intervalListEntry.getKey());
                        
                        SnpIntervalList snpIntervalList = new SnpIntervalList(
                                intervalListEntry.getValue(),
                                newIbsRegions.getStartInBasePairs(),
                                newIbsRegions.getExtentInBasePairs());
                        
                        String strainName = intervalListEntry.getKey();
                        File uniqueFileCache = comparisonStrainsFileMap.get(
                                strainName);
                        if(uniqueFileCache == null)
                        {
                            throw new NullPointerException();
                        }
                        
                        snpBlockListsMap.put(
                                strainName,
                                snpIntervalList);
                        
                        ObjectOutputStream objectOutput = new ObjectOutputStream(
                                new BufferedOutputStream(new FileOutputStream(uniqueFileCache)));
                        objectOutput.writeObject(snpIntervalList);
                        objectOutput.flush();
                        objectOutput.close();
                    }
                }
            }
            
            return new SnpIntervalListGroup(
                    snpBlockListsMap);
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "Failed to cache IBS",
                    ex);
            return null;
        }
    }

    /**
     * Get a cache file to be used for the given parameters. This function
     * will not create the file on disk
     * @param genomeName
     *          the genome name
     * @param strain1
     *          the 1st strain
     * @param strain2
     *          the 2nd strain
     * @param chromosomeNumber
     *          the chromosome number
     * @param minimumExtentInSnps
     *          minimum extent in snps
     * @param minimumExtentInBasePairs
     *          minimum extent in base pairs
     * @return
     *          the file
     */
    private File getCacheFile(
            String genomeName,
            String strain1,
            String strain2,
            int chromosomeNumber,
            int minimumExtentInSnps,
            long minimumExtentInBasePairs)
    {
        final String lesserStrain;
        final String greaterStrain;
        if(strain1.compareTo(strain2) >= 0)
        {
            lesserStrain = strain2;
            greaterStrain = strain1;
        }
        else
        {
            lesserStrain = strain1;
            greaterStrain = strain2;
        }
        
        String cacheHashString =
            genomeName + CONCATINATION_STRING +
            lesserStrain + CONCATINATION_STRING +
            greaterStrain + CONCATINATION_STRING +
            chromosomeNumber + CONCATINATION_STRING +
            minimumExtentInSnps + CONCATINATION_STRING +
            minimumExtentInBasePairs;
        
        try
        {
            cacheHashString = URL_CODEC.encode(cacheHashString);
        }
        catch(EncoderException ex)
        {
            LOG.log(Level.WARNING,
                    "File encoding failed... continuing",
                    ex);
        }
        
        File directory = new File(System.getProperty("java.io.tmpdir"));
        return new File(
                directory,
                cacheHashString);
    }
}
