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

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;

/**
 * Uses a scanning algorithm to find maximal haplotype intervals (maximal in
 * the sense that the intervals are allowed to grow as large as possible).
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntervalScanningHaplotypeEstimator implements HaplotypeEstimator
{
    /**
     * @see #getMinimumConsecutiveSnps()
     */
    private final int minimumConsecutiveSnps;
    
    /**
     * @see #getMinimumNumberOfChromosomes()
     */
    private final int minimumNumberOfChromosomes;
    
    /**
     * Constructor
     * @param minimumConsecutiveSnps
     *          see {@link #getMinimumConsecutiveSnps()}
     * @param minimumNumberOfChromosomes
     *          see {@link #getMinimumNumberOfChromosomes()}
     */
    public IntervalScanningHaplotypeEstimator(
            int minimumConsecutiveSnps,
            int minimumNumberOfChromosomes)
    {
        this.minimumConsecutiveSnps = minimumConsecutiveSnps;
        this.minimumNumberOfChromosomes = minimumNumberOfChromosomes;
    }
    
    /**
     * {@inheritDoc}
     */
    public List<PartitionedInterval> estimateHaplotypeBlocks(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        if(sdpInputStream.getSdpCount() != positionInputStream.getSnpCount())
        {
            throw new IllegalArgumentException(
                    "The number of SDPs should match the number of positions: " +
                    sdpInputStream.getSdpCount() + "!=" +
                    positionInputStream.getSnpCount());
        }
        
        List<PartitionedInterval> haplotypeList =
            new ArrayList<PartitionedInterval>();
        
        Map<BitSet, HaplotypeCandidate> haplotypeCandidateMap =
            new HashMap<BitSet, HaplotypeCandidate>();
        Map<BitSet, HaplotypeCandidate> newHaplotypeCandidateMap =
            new HashMap<BitSet, HaplotypeCandidate>();
        
        // iterate through the SNP positions
        int strainCount = sdpInputStream.getSdpStrainNames().length;
        int chromosomeNumber = positionInputStream.getChromosomeNumber();
        long currPositionBp = -1;
        long prevPositionBp = -1;
        int currSnpIndex = 0;
        while(sdpInputStream.hasNextSdp())
        {
            BitSet currSdp = sdpInputStream.getNextSdp();
            BitSet currSdpCompliment = (BitSet)currSdp.clone();
            currSdpCompliment.flip(0, strainCount);
            currPositionBp = positionInputStream.getNextSnpPositionInBasePairs();
            
            // if either SDP is the empty set then skip to the next snp
            BitSet[] snpGroups = new BitSet[] {
                    currSdp,
                    currSdpCompliment};
            
            Iterator<HaplotypeCandidate> haplotypeCandidateIter =
                haplotypeCandidateMap.values().iterator();
            while(haplotypeCandidateIter.hasNext())
            {
                HaplotypeCandidate currHaplotypeCandidate = haplotypeCandidateIter.next();
                boolean terminateCandidateHaplotype = false;
                boolean atLeastOneIntersection = false;
                
                for(BitSet currSnpGroup: snpGroups)
                {
                    BitSet currHaplotypeCandidateStrains =
                        currHaplotypeCandidate.getStrainBitSet();
                    int currHaplotypeCandidateBitCount =
                        currHaplotypeCandidateStrains.cardinality();
                    
                    // create an intersection of the current snp group and the
                    // current candidate
                    BitSet chromosomeIntersection = (BitSet)currSnpGroup.clone();
                    chromosomeIntersection.and(currHaplotypeCandidateStrains);
                    
                    int intersectionBitCount = chromosomeIntersection.cardinality();
                    assert intersectionBitCount <= currSnpGroup.cardinality();
                    assert intersectionBitCount <= currHaplotypeCandidateBitCount;
                    
                    if(intersectionBitCount > 0)
                    {
                        atLeastOneIntersection = true;
                        
                        // see if all of the haplotype chromosomes are
                        // contained in the current SNP group
                        if(currHaplotypeCandidateBitCount != intersectionBitCount)
                        {
                            // it's a partial intersection.
                            // 1st terminate the old potential haplotype
                            // since it's been broken up
                            terminateCandidateHaplotype = true;
                            
                            // 2nd add the new haplotype to the candidate
                            // list if it is big enough
                            if(intersectionBitCount >= this.minimumNumberOfChromosomes)
                            {
                                // this new potential candidate should go back
                                // as far as the current candidate does
                                int startingSnpIndex =
                                    currHaplotypeCandidate.getStartingSnpIndex();
                                
                                // If there are existing candidates with
                                // matching bit sets then we should not
                                // add this one because they must be
                                // longer
                                if(!haplotypeCandidateMap.containsKey(chromosomeIntersection))
                                {
                                    // if there is a new candidate that
                                    // matches, whoever is further back wins
                                    HaplotypeCandidate matchingCandidate =
                                        newHaplotypeCandidateMap.get(chromosomeIntersection);
                                    if(matchingCandidate == null ||
                                       matchingCandidate.getStartingSnpIndex() > startingSnpIndex)
                                    {
                                        HaplotypeCandidate newCandidate =
                                            new HaplotypeCandidate(
                                                    currHaplotypeCandidate.getStartingPositionBp(),
                                                    startingSnpIndex,
                                                    chromosomeIntersection);
                                        newHaplotypeCandidateMap.put(
                                                chromosomeIntersection,
                                                newCandidate);
                                    }
                                }
                            }
                        }
                    }
                }
                
                // there should be at least one intersection
                assert atLeastOneIntersection;
                
                // see if the haplotype candidate was terminated
                if(terminateCandidateHaplotype)
                {
                    // the reason there is no "+1" in the difference is because
                    // the last valid index was the previous one, not this one
                    int indexDifference =
                        currSnpIndex - currHaplotypeCandidate.getStartingSnpIndex();
                    if(indexDifference >= this.minimumConsecutiveSnps)
                    {
                        // the interval is big enough to put
                        // into the permanent haplotype list
                        long startBp = currHaplotypeCandidate.getStartingPositionBp();
                        PartitionedInterval haplotypeBlock = new PartitionedInterval(
                                chromosomeNumber,
                                startBp,
                                1L + prevPositionBp - startBp,
                                currHaplotypeCandidate.getStrainBitSet());
                        haplotypeList.add(haplotypeBlock);
                    }
                    haplotypeCandidateIter.remove();
                }
            }
            
            // Add the new haplotype candidates then clear them
            haplotypeCandidateMap.putAll(newHaplotypeCandidateMap);
            newHaplotypeCandidateMap.clear();
            
            // Add the SNP groups unless they're already added
            for(BitSet currSnpGroup: snpGroups)
            {
                if(currSnpGroup.cardinality() >= this.minimumNumberOfChromosomes &&
                   !haplotypeCandidateMap.containsKey(currSnpGroup))
                {
                    haplotypeCandidateMap.put(
                            currSnpGroup,
                            new HaplotypeCandidate(
                                    currPositionBp,
                                    currSnpIndex,
                                    currSnpGroup));
                }
            }
            
            // increment the SNP counter and set the previous position
            prevPositionBp = currPositionBp;
            currSnpIndex++;
        }
        
        // Now we just need to clean up all of the haplotype candidates
        for(HaplotypeCandidate currHaplotypeCandidate: haplotypeCandidateMap.values())
        {
            int indexDifference =
                currSnpIndex - currHaplotypeCandidate.getStartingSnpIndex();
            if(indexDifference >= this.minimumConsecutiveSnps)
            {
                long startBp = currHaplotypeCandidate.getStartingPositionBp();
                long extentBp = 1L + currPositionBp - startBp;
                haplotypeList.add(new PartitionedInterval(
                        chromosomeNumber,
                        startBp,
                        extentBp,
                        currHaplotypeCandidate.getStrainBitSet()));
            }
        }
        
        return haplotypeList;
    }
    
    /**
     * the minimum interval size in snps that our algorithm will use to
     * estimate haplotypes
     * @return the minimum consecutive SNP count
     */
    public int getMinimumConsecutiveSnps()
    {
        return this.minimumConsecutiveSnps;
    }
    
    /**
     * Getter for the minimum number of chromosomes that must be in a
     * haplotype before we consider it a haplotype
     * @return the minimumNumberOfChromosomes
     */
    public int getMinimumNumberOfChromosomes()
    {
        return this.minimumNumberOfChromosomes;
    }
    
    /**
     * like a "struct" for holding candidate Haplotype information
     * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
     */
    protected static final class HaplotypeCandidate
    {
        private final long startingPositionBp;
        
        private final int startingSnpIndex;
        
        private final BitSet strainBitSet;
        
        /**
         * Constructor
         * @param startingPositionBp
         *          the starting position in base pairs
         * @param startingSnpIndex
         *          the starting index
         * @param strainBitSet
         *          the candidate strains
         */
        public HaplotypeCandidate(
                long startingPositionBp,
                int startingSnpIndex,
                BitSet strainBitSet)
        {
            this.startingPositionBp = startingPositionBp;
            this.startingSnpIndex = startingSnpIndex;
            this.strainBitSet = strainBitSet;
        }
        
        /**
         * Getter for the starting position in base pairs
         * @return the startingPositionBp
         */
        public long getStartingPositionBp()
        {
            return this.startingPositionBp;
        }

        /**
         * Getter for the starting snp index
         * @return the startingSnpIndex
         */
        public int getStartingSnpIndex()
        {
            return this.startingSnpIndex;
        }

        /**
         * Getter for the candidate chromosomes
         * @return the candidateChromosomes
         */
        public BitSet getStrainBitSet()
        {
            return this.strainBitSet;
        }
        
        /**
         * {@inheritDoc}
         */
        @Override
        public int hashCode()
        {
            return this.startingSnpIndex +
                   this.strainBitSet.hashCode();
        }
        
        /**
         * {@inheritDoc}
         */
        @Override
        public boolean equals(Object otherObject)
        {
            if(otherObject instanceof HaplotypeCandidate)
            {
                HaplotypeCandidate otherCandidate = (HaplotypeCandidate)otherObject;
                return this.startingSnpIndex == otherCandidate.startingSnpIndex &&
                       this.strainBitSet.equals(otherCandidate.strainBitSet);
            }
            else
            {
                return false;
            }
        }
    }
}
