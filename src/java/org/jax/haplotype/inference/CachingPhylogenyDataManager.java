/*
 * Copyright (c) 2008 The Jackson Laboratory
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

package org.jax.haplotype.inference;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.StreamDirection;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.haplotype.phylogeny.inference.IntervalScanner;
import org.jax.haplotype.phylogeny.inference.PhylogenyScanner;
import org.jax.util.datastructure.SequenceUtilities;

/**
 * A phylogeny data manager capable of caching results
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CachingPhylogenyDataManager
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            CachingPhylogenyDataManager.class.getName());
    
    private final IntervalScanner intervalScanner =
        new IntervalScanner();
    private final PhylogenyScanner phylogenyScanner =
        new PhylogenyScanner();
    
    private final CachingGenomeDataManager genomeDataManager;

    private static final CachingPhylogenyDataManager instance =
        new CachingPhylogenyDataManager(
                CachingGenomeDataManager.getInstance());

    protected static final String CONCATINATION_STRING = "&";
    
    private final Map<String, File> resultsCacheFileMap =
        new HashMap<String, File>();

    /**
     * Get the singleton instance
     * @return
     *          the singleton instance
     */
    public static CachingPhylogenyDataManager getInstance()
    {
        return CachingPhylogenyDataManager.instance;
    }
    
    /**
     * private constructor (use {@link #getInstance()} instead)
     * @param genomeDataManager
     *          our source for genotype data
     */
    private CachingPhylogenyDataManager(
            CachingGenomeDataManager genomeDataManager)
    {
        this.genomeDataManager = genomeDataManager;
    }

    /**
     * Get the max-k intervals for the given parameters
     * @param genomeName
     *          the name of the genome
     * @param strainNames
     *          the strain names
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the intervals
     * @throws IOException
     *          if the IO fails
     */
    @SuppressWarnings("unchecked")
    public List<IndexedSnpInterval> getIndexedMaxKIntervals(
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
            throws IOException
    {
        if(LOG.isLoggable(Level.INFO))
        {
            LOG.info("inferring max-k intervals");
        }
        
        synchronized(this)
        {
            File cacheFile = this.getCacheFile(
                    "indexed-interval",
                    genomeName,
                    strainNames,
                    chromosomeNumber);
            
            if(cacheFile.exists())
            {
                ObjectInputStream objectInputStream = new ObjectInputStream(
                        new BufferedInputStream(new FileInputStream(cacheFile)));
                
                try
                {
                    List<IndexedSnpInterval> indexedMaxKIntervals =
                        (List<IndexedSnpInterval>)objectInputStream.readObject();
                    objectInputStream.close();
                    return indexedMaxKIntervals;
                }
                catch(ClassNotFoundException ex)
                {
                    LOG.log(Level.SEVERE,
                            "failed to find class",
                            ex);
                    return null;
                }
            }
            else
            {
                GenomeDataSource genomeDataSource =
                    this.genomeDataManager.getGenomeDataMap().get(genomeName);
                ChromosomeDataSource chromosome = genomeDataSource.getChromosomeDataSources().get(
                        chromosomeNumber);
                
                SdpInputStream forwardStream = chromosome.getSdpInputStream(
                        StreamDirection.FORWARD,
                        strainNames);
                SdpInputStream reverseStream = chromosome.getSdpInputStream(
                        StreamDirection.REVERSE,
                        strainNames);
                SdpInputStream uberStream = chromosome.getSdpInputStream(
                        StreamDirection.FORWARD,
                        strainNames);
                
                List<IndexedSnpInterval> indexedMaxKIntervals = this.intervalScanner.maxKScan(
                        forwardStream,
                        reverseStream,
                        uberStream);
                assert SequenceUtilities.isSorted(indexedMaxKIntervals);
                
                cacheFile.createNewFile();
                cacheFile.deleteOnExit();
                
                ObjectOutputStream objectOutputStream = new ObjectOutputStream(
                        new BufferedOutputStream(new FileOutputStream(cacheFile)));
                objectOutputStream.writeObject(indexedMaxKIntervals);
                objectOutputStream.flush();
                objectOutputStream.close();
                
                return indexedMaxKIntervals;
            }
        }
    }

    /**
     * Get the max-k intervals for the given parameters
     * @param genomeName
     *          the name of the genome
     * @param strainNames
     *          the strain names
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the intervals
     * @throws IOException
     *          if the IO fails
     */
    public List<BasePairInterval> getMaxKIntervals(
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
            throws IOException
    {
        if(LOG.isLoggable(Level.INFO))
        {
            LOG.info("inferring max-k intervals");
        }
        
        synchronized(this)
        {
            File cacheFile = this.getCacheFile(
                    "interval",
                    genomeName,
                    strainNames,
                    chromosomeNumber);
            
            if(cacheFile.exists())
            {
                ObjectInputStream objectInputStream = new ObjectInputStream(
                        new BufferedInputStream(new FileInputStream(cacheFile)));
                
                try
                {
                    List<BasePairInterval> maxKIntervals =
                        (List<BasePairInterval>)objectInputStream.readObject();
                    objectInputStream.close();
                    return maxKIntervals;
                }
                catch(ClassNotFoundException ex)
                {
                    LOG.log(Level.SEVERE,
                            "failed to find class",
                            ex);
                    return null;
                }
            }
            else
            {
                GenomeDataSource genomeDataSource =
                    this.genomeDataManager.getGenomeDataMap().get(genomeName);
                ChromosomeDataSource chromosome = genomeDataSource.getChromosomeDataSources().get(
                        chromosomeNumber);
                
                List<IndexedSnpInterval> indexedMaxKIntervals = this.getIndexedMaxKIntervals(
                        genomeName,
                        strainNames,
                        chromosomeNumber);
                List<BasePairInterval> maxKIntervals = this.intervalScanner.toOrderedPhysicalIntervals(
                        indexedMaxKIntervals,
                        chromosome.getSnpPositionInputStream());
                assert SequenceUtilities.isSorted(maxKIntervals);
                
                cacheFile.createNewFile();
                cacheFile.deleteOnExit();
                
                ObjectOutputStream objectOutputStream = new ObjectOutputStream(
                        new BufferedOutputStream(new FileOutputStream(cacheFile)));
                objectOutputStream.writeObject(maxKIntervals);
                objectOutputStream.flush();
                objectOutputStream.close();
                
                return maxKIntervals;
            }
        }
    }

    /**
     * Get phylogenetic trees for the given parameters
     * @param genomeName
     *          the name of the genome
     * @param strainNames
     *          the strain names
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the phylogeny intervals
     * @throws IOException
     *          if we fail on IO
     * @throws NoValidPhylogenyException
     *          if we run into an interval that doesn't correspond to
     *          a valid perfect phylogeny
     */
    @SuppressWarnings("unchecked")
    public List<PhylogenyTreeNode> getPhylogenies(
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
            throws IOException, NoValidPhylogenyException
    {
        File cacheFile = this.getCacheFile(
                "phylogenies",
                genomeName,
                strainNames,
                chromosomeNumber);
        
        if(cacheFile.exists())
        {
            ObjectInputStream objectInputStream = new ObjectInputStream(
                    new BufferedInputStream(new FileInputStream(cacheFile)));
            
            try
            {
                List<PhylogenyTreeNode> phylogenyIntervals =
                    (List<PhylogenyTreeNode>)objectInputStream.readObject();
                objectInputStream.close();
                return phylogenyIntervals;
            }
            catch(ClassNotFoundException ex)
            {
                LOG.log(Level.SEVERE,
                        "failed to find class",
                        ex);
                return null;
            }
        }
        else
        {
            GenomeDataSource genomeDataSource =
                this.genomeDataManager.getGenomeDataMap().get(genomeName);
            ChromosomeDataSource chromosome = genomeDataSource.getChromosomeDataSources().get(
                    chromosomeNumber);
            
            List<IndexedSnpInterval> indexedMaxKIntervals = this.getIndexedMaxKIntervals(
                    genomeName,
                    strainNames,
                    chromosomeNumber);
            List<PhylogenyTreeNode> phylogenies = this.phylogenyScanner.inferPerfectPhylogenies(
                    chromosome.getSdpInputStream(strainNames),
                    indexedMaxKIntervals);
            assert indexedMaxKIntervals.size() == phylogenies.size();
            
            cacheFile.createNewFile();
            cacheFile.deleteOnExit();
            
            ObjectOutputStream objectOutputStream = new ObjectOutputStream(
                    new BufferedOutputStream(new FileOutputStream(cacheFile)));
            objectOutputStream.writeObject(phylogenies);
            objectOutputStream.flush();
            objectOutputStream.close();
            
            return phylogenies;
        }
    }
    
    /**
     * Get phylogenetic trees for the given parameters
     * @param genomeName
     *          the name of the genome
     * @param strainNames
     *          the strain names
     * @param chromosomeNumber
     *          the chromosome number
     * @return
     *          the phylogeny intervals
     * @throws IOException
     *          if we fail on IO
     * @throws NoValidPhylogenyException
     *          if we run into an interval that doesn't correspond to
     *          a valid perfect phylogeny
     */
    @SuppressWarnings("unchecked")
    public List<PhylogenyInterval> getPhylogeneticIntervals(
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
            throws IOException, NoValidPhylogenyException
    {
        File cacheFile = this.getCacheFile(
                "phylo-intervals",
                genomeName,
                strainNames,
                chromosomeNumber);
        
        if(cacheFile.exists())
        {
            ObjectInputStream objectInputStream = new ObjectInputStream(
                    new BufferedInputStream(new FileInputStream(cacheFile)));
            
            try
            {
                List<PhylogenyInterval> phylogenyIntervals =
                    (List<PhylogenyInterval>)objectInputStream.readObject();
                objectInputStream.close();
                return phylogenyIntervals;
            }
            catch(ClassNotFoundException ex)
            {
                LOG.log(Level.SEVERE,
                        "failed to find class",
                        ex);
                return null;
            }
        }
        else
        {
            GenomeDataSource genomeDataSource =
                this.genomeDataManager.getGenomeDataMap().get(genomeName);
            ChromosomeDataSource chromosome = genomeDataSource.getChromosomeDataSources().get(
                    chromosomeNumber);
            
            List<IndexedSnpInterval> indexedMaxKIntervals = this.getIndexedMaxKIntervals(
                    genomeName,
                    strainNames,
                    chromosomeNumber);
            List<PhylogenyTreeNode> phylogenies = this.phylogenyScanner.inferPerfectPhylogenies(
                    chromosome.getSdpInputStream(strainNames),
                    indexedMaxKIntervals);
            List<BasePairInterval> maxKIntervals = this.intervalScanner.toOrderedPhysicalIntervals(
                    indexedMaxKIntervals,
                    chromosome.getSnpPositionInputStream());
            int phylogenyCount = phylogenies.size();
            assert maxKIntervals.size() == phylogenyCount;
            
            List<PhylogenyInterval> phylogenyIntervals = new ArrayList<PhylogenyInterval>(
                    phylogenyCount);
            for(int i = 0; i < phylogenyCount; i++)
            {
                phylogenyIntervals.add(new PhylogenyInterval(
                        phylogenies.get(i),
                        maxKIntervals.get(i)));
            }
            
            cacheFile.createNewFile();
            cacheFile.deleteOnExit();
            
            ObjectOutputStream objectOutputStream = new ObjectOutputStream(
                    new BufferedOutputStream(new FileOutputStream(cacheFile)));
            objectOutputStream.writeObject(phylogenyIntervals);
            objectOutputStream.flush();
            objectOutputStream.close();
            
            return phylogenyIntervals;
        }
    }
    
    /**
     * Get the genome data manager
     * @return
     *          the genome data manager
     */
    public CachingGenomeDataManager getGenomeDataManager()
    {
        return this.genomeDataManager;
    }

    /**
     * Get a cache file to be used for the given parameters. This function
     * will not create the file on disk
     */
    private File getCacheFile(
            String filePrefix,
            String genomeName,
            String[] strainNames,
            int chromosomeNumber)
    {
        strainNames = strainNames.clone();
        Arrays.sort(strainNames);
        
        String cacheKeyString =
            filePrefix + CONCATINATION_STRING +
            genomeName + CONCATINATION_STRING +
            Arrays.toString(strainNames) + CONCATINATION_STRING +
            chromosomeNumber;
        
        File cacheFile = this.resultsCacheFileMap.get(cacheKeyString);
        if(cacheFile == null)
        {
            File directory = new File(System.getProperty("java.io.tmpdir"));
            cacheFile = new File(
                    directory,
                    filePrefix + "-" + this.resultsCacheFileMap.size());
            this.resultsCacheFileMap.put(
                    cacheKeyString,
                    cacheFile);
        }
        return cacheFile;
    }
}
