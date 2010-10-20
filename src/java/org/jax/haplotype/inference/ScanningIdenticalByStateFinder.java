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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.SnpIntervalList;
import org.jax.geneticutil.data.SnpIntervalListGroup;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.util.datastructure.SequenceUtilities;

/**
 * A scanning (front to back) implementation for finding chromosome
 * IBS regions
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ScanningIdenticalByStateFinder implements IdenticalByStateFinder
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            ScanningIdenticalByStateFinder.class.getName());
    
    /**
     * {@inheritDoc}
     */
    public SnpIntervalListGroup findIdenticalByStateRegions(
            final StrainChromosome referenceStrainChromosome,
            final Set<StrainChromosome> comparisonStrainChromosomes,
            final int minimumExtentInSnps,
            final long minimumExtentInBasePairs)
    {
        Map<String, List<BasePairInterval>> ibsBlocksMap =
            new HashMap<String, List<BasePairInterval>>();
        
        long minStartPositionInBasePairs = Long.MAX_VALUE;
        long maxEndPositionInBasePairs = 0L;
        for(StrainChromosome chromosome: comparisonStrainChromosomes)
        {
            SnpIntervalList ibsList = this.findIdenticalByStateRegions(
                    referenceStrainChromosome,
                    chromosome,
                    minimumExtentInSnps,
                    minimumExtentInBasePairs);
            long newMinPositionInBasePairs = ibsList.getStartInBasePairs();
            long newEndPositionInBasePairs =
                newMinPositionInBasePairs +
                ibsList.getExtentInBasePairs();
            
            if(newMinPositionInBasePairs < minStartPositionInBasePairs)
            {
                minStartPositionInBasePairs = newMinPositionInBasePairs;
            }
            
            if(newEndPositionInBasePairs > maxEndPositionInBasePairs)
            {
                maxEndPositionInBasePairs = newEndPositionInBasePairs;
            }
            
            ibsBlocksMap.put(
                    chromosome.getStrainName(),
                    ibsList.getSnpBlocks());
        }
        
        return new SnpIntervalListGroup(
                ibsBlocksMap,
                minStartPositionInBasePairs,
                maxEndPositionInBasePairs - minStartPositionInBasePairs);
    }
    
    /**
     * find ibs regions using streaming data types
     * @param comparisonSdpStream
     *          the comparison strain data
     * @param snpPositionInputStream
     *          the SNP position data
     * @param minimumExtentInSnps
     *          the minimum number of consecutive SNPs before we call something
     *          IBS
     * @param minimumExtentInBasePairs
     *          the minimum extent in base pairs before we call something IBS
     * @return
     *          the intervals
     */
    public SnpIntervalListGroup findIdenticalByStateRegions(
            final SdpInputStream comparisonSdpStream,
            final SnpPositionInputStream snpPositionInputStream,
            final long minimumExtentInSnps,
            final long minimumExtentInBasePairs)
    {
        try
        {
            final String[] comparisonStrainNames = comparisonSdpStream.getSdpStrainNames();
            
            long[] lastMatchingSnpIndices = new long[comparisonStrainNames.length];
            Arrays.fill(lastMatchingSnpIndices, -1L);
            long[] lastMatchingSnpPositions = new long[comparisonStrainNames.length];
            Arrays.fill(lastMatchingSnpPositions, -1L);
            
            List<BasePairInterval>[] snpIntervalLists =
                SequenceUtilities.instantiateGenericArray(
                        List.class,
                        comparisonStrainNames.length);
            for(int i = 0; i < snpIntervalLists.length; i++)
            {
                snpIntervalLists[i] = new ArrayList<BasePairInterval>();
            }
            
            long snpIndex = 0L;
            long prevReferencePosition = -1L;
            long nextReferencePosition;
            while(comparisonSdpStream.hasNextSdp())
            {
                nextReferencePosition =
                    snpPositionInputStream.getNextSnpPositionInBasePairs();
                BitSet nextSdp = comparisonSdpStream.getNextSdp();
                for(int i = 0; i < comparisonStrainNames.length; i++)
                {
                    long lastMatchingSnpIndex = lastMatchingSnpIndices[i];
                    if(nextSdp.get(i))
                    {
                        if(lastMatchingSnpIndex == -1L)
                        {
                            lastMatchingSnpIndices[i] = snpIndex;
                            lastMatchingSnpPositions[i] = nextReferencePosition;
                        }
                    }
                    else
                    {
                        if(lastMatchingSnpIndex >= 0)
                        {
                            long lastMatchingSnpPosition =
                                lastMatchingSnpPositions[i];
                            long extentInSnps =
                                snpIndex - lastMatchingSnpIndex;
                            long extentInBasePairs =
                                1L + prevReferencePosition - lastMatchingSnpPosition;
                            if(extentInBasePairs >= minimumExtentInBasePairs &&
                               extentInSnps >= minimumExtentInSnps)
                            {
                                snpIntervalLists[i].add(new SimpleBasePairInterval(
                                        snpPositionInputStream.getChromosomeNumber(),
                                        lastMatchingSnpPosition,
                                        extentInBasePairs));
                            }
                        }
                        
                        lastMatchingSnpIndices[i] = -1L;
                        lastMatchingSnpPositions[i] = -1L;
                    }
                }
                
                snpIndex++;
                prevReferencePosition = nextReferencePosition;
            }
            
            // cleanup
            for(int i = 0; i < comparisonStrainNames.length; i++)
            {
                long lastMatchingSnpIndex = lastMatchingSnpIndices[i];
                if(lastMatchingSnpIndex >= 0)
                {
                    long lastMatchingSnpPosition =
                        lastMatchingSnpPositions[i];
                    long extentInSnps =
                        snpIndex - lastMatchingSnpIndex;
                    long extentInBasePairs =
                        1L + prevReferencePosition - lastMatchingSnpPosition;
                    if(extentInBasePairs >= minimumExtentInBasePairs &&
                       extentInSnps >= minimumExtentInSnps)
                    {
                        snpIntervalLists[i].add(new SimpleBasePairInterval(
                                snpPositionInputStream.getChromosomeNumber(),
                                lastMatchingSnpPosition,
                                extentInBasePairs));
                    }
                }
            }
            
            Map<String, List<BasePairInterval>> intervalMap =
                new HashMap<String, List<BasePairInterval>>(comparisonStrainNames.length);
            for(int i = 0; i < comparisonStrainNames.length; i++)
            {
                intervalMap.put(
                        comparisonStrainNames[i],
                        snpIntervalLists[i]);
            }
            
            SnpIntervalListGroup snpIntervalListGroup =
                new SnpIntervalListGroup(
                        intervalMap,
                        snpPositionInputStream.getStartInBasePairs(),
                        snpPositionInputStream.getExtentInBasePairs());
            return snpIntervalListGroup;
        }
        catch(IOException ex)
        {
            LOG.log(Level.SEVERE, "failed to find IBS regions", ex);
            return null;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    public SnpIntervalList findIdenticalByStateRegions(
            final StrainChromosome chromosome1,
            final StrainChromosome chromosome2,
            final int minimumExtentInSnps,
            final long minimumExtentInBasePairs)
    {
        System.out.println("calling scanning findIdenticalByStateRegions");
        SingleNucleotidePolymorphism[] snps1 =
            chromosome1.getSingleNucleotidePolymorphisms();
        SingleNucleotidePolymorphism[] snps2 =
            chromosome2.getSingleNucleotidePolymorphisms();
        
        if(chromosome1.getChromosomeNumber() != chromosome2.getChromosomeNumber())
        {
            throw new IllegalArgumentException(
                    "cant find IBS for different chromosome numbers");
        }
        else if(snps1.length != snps2.length)
        {
            throw new IllegalArgumentException(
                    "IBS algorithm requires that # of SNPS is identical" +
                    "between chromosomes");
        }
        
        List<BasePairInterval> ibsBlocks =
            new ArrayList<BasePairInterval>();
        int matchingIndexStart = -1;
        for(int snpIndex = 0; snpIndex < snps1.length; snpIndex++)
        {
            SingleNucleotidePolymorphism snp1 = snps1[snpIndex];
            SingleNucleotidePolymorphism snp2 = snps2[snpIndex];
            
            if(snp1.getSnpType() == snp2.getSnpType())
            {
                if(matchingIndexStart == -1)
                {
                    matchingIndexStart = snpIndex;
                }
            }
            else
            {
                if(matchingIndexStart != -1)
                {
                    int snpExtent = snpIndex - matchingIndexStart;
                    long basePairExtent =
                        1 +
                        (snps1[snpIndex - 1].getPositionInBasePairs() -
                        snps1[matchingIndexStart].getPositionInBasePairs());
                    if(snpExtent >= minimumExtentInSnps &&
                       basePairExtent >= minimumExtentInBasePairs)
                    {
                        BasePairInterval block =
                            new SimpleBasePairInterval(
                                    chromosome1.getChromosomeNumber(),
                                    snps1[matchingIndexStart].getPositionInBasePairs(),
                                    basePairExtent);
                        ibsBlocks.add(block);
                    }
                    
                    matchingIndexStart = -1;
                }
            }
        }
        
        // clean up
        if(matchingIndexStart != -1)
        {
            int snpExtent = snps1.length - matchingIndexStart;
            long basePairExtent =
                (snps1[snps1.length - 1].getPositionInBasePairs() -
                snps1[matchingIndexStart].getPositionInBasePairs()) +
                1;
            if(snpExtent >= minimumExtentInSnps &&
               basePairExtent >= minimumExtentInBasePairs)
            {
                BasePairInterval block =
                    new SimpleBasePairInterval(
                            chromosome1.getChromosomeNumber(),
                            snps1[matchingIndexStart].getPositionInBasePairs(),
                            basePairExtent);
                ibsBlocks.add(block);
            }
        }
        
        return new SnpIntervalList(
                ibsBlocks,
                snps1[0].getPositionInBasePairs(),
                snps1[snps1.length - 1].getPositionInBasePairs() -
                snps1[0].getPositionInBasePairs() +
                1L);
    }
}
