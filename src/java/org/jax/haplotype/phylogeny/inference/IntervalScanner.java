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

package org.jax.haplotype.phylogeny.inference;

import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.haplotype.io.MinorityNormalizedSdpInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.haplotype.io.StreamDirection;
import org.jax.util.datastructure.ExposedArrayList;
import org.jax.util.datastructure.SequenceUtilities;

/**
 * A class for building compatible intervals where a compatible interval
 * is a contiguous region of the genome where all SNP value differences can
 * be explained by mutations using the "infinite sites" assumption. One
 * nicety of compatible intervals is that you can always build perfect
 * phylogenies from them.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntervalScanner
{
    /**
     * Class for grouping an SDP with its index
     */
    private static class SdpIndexPair
    {
        private final BitSet sdpBits;
        
        private final int index;

        /**
         * Constructor
         * @param sdpBits
         *          see {@link #getSdpBits()}
         * @param index
         *          see {@link #getIndex()}
         */
        public SdpIndexPair(BitSet sdpBits, int index)
        {
            this.sdpBits = sdpBits;
            this.index = index;
        }
        
        /**
         * Getter for the index
         * @return the index
         */
        public int getIndex()
        {
            return this.index;
        }
        
        /**
         * getter for the SDP
         * @return the SDP
         */
        public BitSet getSdpBits()
        {
            return this.sdpBits;
        }
    }
    
    /**
     * Do a max-k scan which involves doing a lot of other scans.
     * All of the streams passed in should basically represent the same
     * underlying data (except the reverse stream should read in the
     * {@link StreamDirection#REVERSE} direction
     * @param forwardStream
     *          the stream that we do a greedy forward scan on
     * @param reverseStream
     *          the stream that we do a reverse scan on
     * @param uberStream
     *          the stream that we do an uber scan on
     * @return
     *          the max-k interval
     * @throws IOException
     *          if the streams throw an exception
     */
    public List<IndexedSnpInterval> maxKScan(
            SdpInputStream forwardStream,
            SdpInputStream reverseStream,
            SdpInputStream uberStream) throws IOException
    {
        if(forwardStream.getReadDirection() != StreamDirection.FORWARD)
        {
            throw new IllegalArgumentException(
                    "forwardStream must read forward");
        }
        
        if(reverseStream.getReadDirection() != StreamDirection.REVERSE)
        {
            throw new IllegalArgumentException(
                    "reverseStream must read in reverse");
        }
        
        List<IndexedSnpInterval> forwardIntervals =
            this.greedyScan(forwardStream);
        List<IndexedSnpInterval> reverseIntervals =
            this.greedyScan(reverseStream);
        List<IndexedSnpInterval> uberIntervals =
            this.uberScan(uberStream);
        List<IndexedSnpInterval> coreIntervals =
            this.createCoreIntervals(forwardIntervals, reverseIntervals);
        List<List<IndexedSnpInterval>> uberCores =
            this.createUberCores(uberIntervals, coreIntervals);
        List<IndexedSnpInterval> maxKIntervals =
            this.createMaxKIntervals(uberCores);
        
        return maxKIntervals;
    }
    
    /**
     * Do an uber-scan looking for every possible maximal compatible interval
     * @param sdpInputStream
     *          an SDP input stream to read from
     * @return
     *          the uber list of compatible intervals
     * @throws IOException
     *          if we get an exception reading from the SDP stream
     */
    public List<IndexedSnpInterval> uberScan(
            SdpInputStream sdpInputStream)
            throws IOException
    {
        if(sdpInputStream.getReadDirection() == StreamDirection.REVERSE)
        {
            throw new IllegalArgumentException(
                    "uber scan only works on forward streams");
        }
        
        // this algorithm assumes minority normalized SDPs
        sdpInputStream = new MinorityNormalizedSdpInputStream(sdpInputStream);
        
        ArrayList<IndexedSnpInterval> intervals = new ArrayList<IndexedSnpInterval>();
        
        int currSdpIndex = -1;
        int startIndex = 0;
        ExposedArrayList<SdpIndexPair> intervalSdps =
            new ExposedArrayList<SdpIndexPair>();
        while(sdpInputStream.hasNextSdp())
        {
            BitSet currSdp = sdpInputStream.getNextSdp();
            currSdpIndex++;
            
            int nearestIncompatibleIndex = -1;
            while(sdpInputStream.hasNextSdp() && nearestIncompatibleIndex == -1)
            {
                nearestIncompatibleIndex = this.testCompatibleAndUberAdd(
                        intervalSdps,
                        currSdp,
                        currSdpIndex);
                
                if(nearestIncompatibleIndex == -1)
                {
                    currSdp = sdpInputStream.getNextSdp();
                    currSdpIndex++;
                }
            }
            
            if(sdpInputStream.hasNextSdp())
            {
                // we know that we hit an
                // incompatibility so we add the interval
                assert nearestIncompatibleIndex >= 0;
                intervals.add(new IndexedSnpInterval(
                        startIndex,
                        currSdpIndex - startIndex));
                
                // reinitialize the start index position (we start just after
                // the nearest incompatible SDP
                startIndex = nearestIncompatibleIndex + 1;
            }
            else
            {
                // we're done streaming. clean up by adding final interval[s]
                nearestIncompatibleIndex = this.testCompatibleAndUberAdd(
                        intervalSdps,
                        currSdp,
                        currSdpIndex);
                if(nearestIncompatibleIndex >= 0)
                {
                    // The final SDP is incompatible, so
                    // add the last 2 intervals
                    intervals.add(new IndexedSnpInterval(
                            startIndex,
                            currSdpIndex - startIndex));
                    
                    int furthestCompatibleIndex = nearestIncompatibleIndex + 1;
                    assert furthestCompatibleIndex <= currSdpIndex;
                    intervals.add(new IndexedSnpInterval(
                            furthestCompatibleIndex,
                            1 + currSdpIndex - furthestCompatibleIndex));
                }
                else
                {
                    // The final SDP is compatible, so
                    // add the last interval
                    intervals.add(new IndexedSnpInterval(
                            startIndex,
                            1 + currSdpIndex - startIndex));
                }
            }
        }
        
        intervals.trimToSize();
        
        return intervals;
    }
    
    /**
     * Create a set of core intervals from the given forward and reverse
     * greedy scan results
     * @param forwardGreedyIntervals
     *          the forward scan intervals
     * @param reverseGreedyIntervals
     *          the reverse scan intervals
     * @return
     *          the cores
     */
    public List<IndexedSnpInterval> createCoreIntervals(
            List<IndexedSnpInterval> forwardGreedyIntervals,
            List<IndexedSnpInterval> reverseGreedyIntervals)
    {
        if(forwardGreedyIntervals.size() != reverseGreedyIntervals.size())
        {
            throw new IllegalArgumentException(
                    "the reverse and forward interval lists should be the " +
                    "same size");
        }
        
        int intervalCount = forwardGreedyIntervals.size();
        List<IndexedSnpInterval> coreIntervals =
            new ArrayList<IndexedSnpInterval>(intervalCount);
        for(int intervalIndex = 0; intervalIndex < intervalCount; intervalIndex++)
        {
            IndexedSnpInterval currForwardInterval =
                forwardGreedyIntervals.get(intervalIndex);
            IndexedSnpInterval currReverseInterval =
                reverseGreedyIntervals.get(intervalIndex);
            
            int coreIntervalStart = currForwardInterval.getStartIndex();
            int coreIntervalEnd = currReverseInterval.getEndIndex();
            assert coreIntervalStart <= coreIntervalEnd;
            
            coreIntervals.add(new IndexedSnpInterval(
                    coreIntervalStart,
                    1 + coreIntervalEnd - coreIntervalStart));
        }
        
        return coreIntervals;
    }
    
    /**
     * Prune the uber intervals in preparation for calculating the max-k
     * intervals. We can do this by throwing out any uber intervals that do
     * not intersect exactly one core interval since we know that max-k
     * intervals should intersect one and only one core interval
     * @param uberIntervals
     *          the uber set of intervals (unmodified by this function)
     * @param coreIntervals
     *          the core set of intervals (unmodified by this function)
     * @return
     *          a new list of intervals which is a subset of the uber set
     *          following some rules that we know apply to the max-k set
     *          of intervals
     */
    public List<List<IndexedSnpInterval>> createUberCores(
            List<IndexedSnpInterval> uberIntervals,
            List<IndexedSnpInterval> coreIntervals)
    {
        if(uberIntervals.size() < coreIntervals.size())
        {
            throw new IllegalArgumentException(
                    "the list of uber intervals should be at least as big as " +
                    "the list of core intervals");
        }
        
        int uberSize = uberIntervals.size();
        int coreSize = coreIntervals.size();
        List<List<IndexedSnpInterval>> uberIntervalsSubset =
            new ArrayList<List<IndexedSnpInterval>>(coreSize);
        
        if(coreSize >= 1)
        {
            // initialize the data that we'll use in our subsetting loop
            int coreIndex = 0;
            IndexedSnpInterval prevCore = null;
            IndexedSnpInterval currCore = coreIntervals.get(coreIndex);
            IndexedSnpInterval nextCore = null;
            if(coreIndex + 1 < coreSize)
            {
                nextCore = coreIntervals.get(coreIndex + 1);
            }
            
            ArrayList<IndexedSnpInterval> currCoreUberIntervals =
                new ArrayList<IndexedSnpInterval>();
            
            // OK, take care of the subsetting
            for(int uberIndex = 0; uberIndex < uberSize && currCore != null; uberIndex++)
            {
                IndexedSnpInterval currUberInterval =
                    uberIntervals.get(uberIndex);
                
                int uberStart = currUberInterval.getStartIndex();
                
                // do we need to move on to the next core?
                if(uberStart > currCore.getEndIndex())
                {
                    assert !currCoreUberIntervals.isEmpty();
                    assert SequenceUtilities.isSorted(currCoreUberIntervals);
                    
                    // move on to the next core
                    coreIndex++;
                    prevCore = currCore;
                    currCore = nextCore;
                    nextCore = null;
                    if(coreIndex + 1 < coreSize)
                    {
                        nextCore = coreIntervals.get(coreIndex + 1);
                    }
                    uberIntervalsSubset.add(currCoreUberIntervals);
                    currCoreUberIntervals.trimToSize();
                    currCoreUberIntervals =
                        new ArrayList<IndexedSnpInterval>(1);
                }
                
                if(currCore != null)
                {
                    assert uberStart <= currCore.getEndIndex();
                    
                    // we have the correct core, now see if the uber interval
                    // has any chance to become a max-k interval (for this it
                    // must cover only the current core, without intersecting
                    // with the previous or next cores)
                    if(currUberInterval.contains(currCore) &&
                       (prevCore == null || !currUberInterval.intersects(prevCore)) &&
                       (nextCore == null || !currUberInterval.intersects(nextCore)))
                    {
                        // pass this one through
                        currCoreUberIntervals.add(currUberInterval);
                    }
                }
            }
            
            // clean-up
            if(!currCoreUberIntervals.isEmpty())
            {
                assert SequenceUtilities.isSorted(currCoreUberIntervals);
                uberIntervalsSubset.add(currCoreUberIntervals);
            }
            
            assert uberIntervalsSubset.size() == coreSize;
        }
        
        return uberIntervalsSubset;
    }

    /**
     * Add to the interval SDPs assuming that we're working toward the
     * Uber set of compatible intervals. This function modifies the list
     * in different ways depending on whether it finds a match or a conflict,
     * but you shouldn't really have to care about that.
     * <br/><br/>
     * If this function finds a conflict, it returns the the index of the
     * nearest conflicting SNP.
     * @param intervalSdps
     *          the list of SDPs for the current interval. this list is "owned"
     *          by this function and it's probably best not to touch it or try
     *          to make sense of it outside of this function
     * @param sdpToAdd
     *          the SDP that we're trying to introduce to the interval
     * @param sdpIndex
     *          the index of the SDP
     * @return
     *          -1 if the given SDP is compatible with current interval SDPs
     *          or the index of the nearest incompatibility if we find one
     */
    private int testCompatibleAndUberAdd(
            ExposedArrayList<SdpIndexPair> intervalSdps,
            BitSet sdpToAdd,
            int sdpIndex)
    {
        int intervalSdpCount = intervalSdps.size();
        for(int i = intervalSdpCount - 1; i >= 0 ; i--)
        {
            SdpIndexPair currPair = intervalSdps.get(i);
            BitSet currSdp = currPair.getSdpBits();
            if(sdpToAdd.equals(currSdp))
            {
                // this SDP was already added, so return compatible after
                // moving the snp to the end of the list
                intervalSdps.remove(i);
                intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
                return -1;
            }
            else if(!areMinorityNormalizedSdpsCompatible(sdpToAdd, currSdp))
            {
                // found an incompatibility. remove everything before the
                // incompatibility and return the index of the sdp that
                // caused the trouble
                intervalSdps.removeRange(
                        0,
                        i + 1);
                intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
                
                return currPair.getIndex();
            }
        }
        
        // we're compatible and didn't find any exact matches for the SDP
        // which means that we need to add it to the list
        intervalSdps.add(new SdpIndexPair(sdpToAdd, sdpIndex));
        return -1;
    }
    
    /**
     * Perform a greedy scan on the given streams. If the stream we're given
     * reads in the reverse direction then we'll do a
     * {@link #reverseIndexedIntervals(List, int)} on the intervals before
     * returning. This means that the indices should be easily comparable with
     * intervals read in the forward direction.
     * @param sdpInputStream
     *          the SDP's to scan
     * @return
     *          the list of SNP intervals. these intervals will cover the
     *          entire genome with no overlap
     * @throws IOException
     *          if we run into unexpected trouble reading the stream
     */
    public List<IndexedSnpInterval> greedyScan(
            SdpInputStream sdpInputStream)
            throws IOException
    {
        // this algorithm assumes minority normalized SDPs
        sdpInputStream = new MinorityNormalizedSdpInputStream(sdpInputStream);
        
        ArrayList<IndexedSnpInterval> intervals = new ArrayList<IndexedSnpInterval>();
        
        int currSdpIndex = -1;
        int startIndex = 0;
        List<BitSet> intervalSdps = new ArrayList<BitSet>();
        while(sdpInputStream.hasNextSdp())
        {
            BitSet currSdp = sdpInputStream.getNextSdp();
            currSdpIndex++;
            while(sdpInputStream.hasNextSdp() &&
                  this.checkCompatibilityAndAddSdp(intervalSdps, currSdp))
            {
                currSdp = sdpInputStream.getNextSdp();
                currSdpIndex++;
            }
            
            if(sdpInputStream.hasNextSdp())
            {
                // since hasNextSdp is true we know that we hit an
                // incompatibility so we add the interval
                intervals.add(new IndexedSnpInterval(
                        startIndex,
                        currSdpIndex - startIndex));
                
                // reinitialize for next interval (we need to include conflict SDP)
                intervalSdps.clear();
                intervalSdps.add(currSdp);
                startIndex = currSdpIndex;
            }
            else
            {
                // cleanup by adding final interval[s]
                if(!this.checkCompatibilityAndAddSdp(intervalSdps, currSdp))
                {
                    // The final SDP is incompatible, so
                    // add the last 2 intervals
                    intervals.add(new IndexedSnpInterval(
                            startIndex,
                            currSdpIndex - startIndex));
                    intervals.add(new IndexedSnpInterval(
                            currSdpIndex,
                            1));
                }
                else
                {
                    // The final SDP is compatible, so
                    // add the last interval
                    intervals.add(new IndexedSnpInterval(
                            startIndex,
                            1 + currSdpIndex - startIndex));
                }
            }
        }
        
        intervals.trimToSize();
        
        // we promised to reverse the intervals if we're reading a reverse
        // stream
        if(sdpInputStream.getReadDirection() == StreamDirection.REVERSE)
        {
            assert sdpInputStream.getSdpCount() < Integer.MAX_VALUE;
            this.reverseIndexedIntervals(
                    intervals,
                    (int)sdpInputStream.getSdpCount());
        }
        
        return intervals;
    }

    /**
     * Check the compatibility of the SDP and if it's compatible, add it
     * @param intervalSdps
     *          the interval SDPs to compare against
     * @param sdpToAdd
     *          the SDP we'll try to add to the interval
     * @return
     *          true if the SDP is fully compatible, false otherwise
     */
    private boolean checkCompatibilityAndAddSdp(
            List<BitSet> intervalSdps,
            BitSet sdpToAdd)
    {
        int intervalSdpCount = intervalSdps.size();
        for(int i = 0; i < intervalSdpCount; i++)
        {
            BitSet currSdp = intervalSdps.get(i);
            if(sdpToAdd.equals(currSdp))
            {
                // this SDP was already added, so return compatible
                // without adding SDP
                return true;
            }
            else if(!areMinorityNormalizedSdpsCompatible(sdpToAdd, currSdp))
            {
                // found an incompatibility, so return incompatible
                // without adding SDP
                return false;
            }
        }
        
        // we're compatible and didn't find any exact matches for the SDP
        // which means that we need to add it to the list
        intervalSdps.add(sdpToAdd);
        return true;
    }
    
    /**
     * Test if the two SDPs are compatible assuming that they follow the
     * normalization rules from {@link MinorityNormalizedSdpInputStream}
     * @param minorityNormalizedSdp1
     *          the 1st SDP
     * @param minorityNormalizedSdp2
     *          the 2nd SDP
     * @return
     *          true if the given SDPs are compatible
     */
    public static boolean areMinorityNormalizedSdpsCompatible(
            BitSet minorityNormalizedSdp1,
            BitSet minorityNormalizedSdp2)
    {
        if(!minorityNormalizedSdp1.intersects(minorityNormalizedSdp2))
        {
            // disjoint SDPs are compatible
            return true;
        }
        else
        {
            BitSet intersection = (BitSet)minorityNormalizedSdp1.clone();
            intersection.and(minorityNormalizedSdp2);
            
            if(intersection.equals(minorityNormalizedSdp1))
            {
                // SDP1 is a subset of SDP2 indicating compatible
                return true;
            }
            else if(intersection.equals(minorityNormalizedSdp2))
            {
                // SDP2 is a subset of SDP1 indicating compatible
                return true;
            }
            else
            {
                // they're incompatible since they're not disjoint and they
                // don't have a subset/superset relationship
                return false;
            }
        }
    }
    
    /**
     * Performs {@link #reverseIndexedInterval(IndexedSnpInterval, int)} on
     * every interval in the list and also reverses the ordering of intervals
     * in the list
     * @param intervalsToReverse
     *          the list to reverse
     * @param totalSnpCount
     *          the total snp count (range) that we're flipping the intervals
     *          over
     */
    protected void reverseIndexedIntervals(
            List<IndexedSnpInterval> intervalsToReverse,
            int totalSnpCount)
    {
        int intervalCount = intervalsToReverse.size();
        int halfIntervalCount = intervalCount / 2;
        
        for(int index1 = 0; index1 < halfIntervalCount; index1++)
        {
            int index2 = (intervalCount - index1) - 1;
            
            IndexedSnpInterval interval1 = intervalsToReverse.get(index1);
            IndexedSnpInterval interval2 = intervalsToReverse.get(index2);
            
            intervalsToReverse.set(
                    index1,
                    this.reverseIndexedInterval(interval2, totalSnpCount));
            intervalsToReverse.set(
                    index2,
                    this.reverseIndexedInterval(interval1, totalSnpCount));
        }
        
        // if we have an odd count we still need to reverse the middle interval
        if(intervalCount % 2 == 1)
        {
            int middleIndex = halfIntervalCount;
            IndexedSnpInterval middleInterval =
                intervalsToReverse.get(middleIndex);
            intervalsToReverse.set(
                    middleIndex,
                    this.reverseIndexedInterval(middleInterval, totalSnpCount));
        }
    }
    
    /**
     * Reverse the indecies assuming the given total index count
     * @param intervalToReverse
     *          the interval 
     * @param totalSnpCount
     *          the total # of SNPs. we need to know this in order to be able
     *          to flip the interval
     * @return
     *          the interval which will be reveresed (or "flipped") along the
     *          range indicated by the total count
     */
    protected IndexedSnpInterval reverseIndexedInterval(
            IndexedSnpInterval intervalToReverse,
            int totalSnpCount)
    {
        int extent = intervalToReverse.getExtentInIndices();
        int newStart =
            (totalSnpCount - intervalToReverse.getStartIndex()) - extent;
        return new IndexedSnpInterval(
                newStart,
                extent);
    }
    
    /**
     * Scan the given uber cores to come up with the max-k data set
     * @param uberCores
     *          the uber cores to search. see
     *          {@link #createUberCores(List, List)}
     * @return
     *          the max-k interval list
     */
    public List<IndexedSnpInterval> createMaxKIntervals(List<List<IndexedSnpInterval>> uberCores)
    {
        int coreCount = uberCores.size();
        List<IndexedSnpInterval> maxKIntervals = new ArrayList<IndexedSnpInterval>(
                coreCount);
        
        if(coreCount >= 1)
        {
            // initialize data pre-scan
            int[][] forwardPointers = new int[uberCores.size() - 1][];
            List<IndexedSnpInterval> prevCoreGroup = uberCores.get(
                    uberCores.size() - 1);
            long[] cumulativeExtents = new long[prevCoreGroup.size()];
            for(int i = 0; i < cumulativeExtents.length; i++)
            {
                cumulativeExtents[i] = prevCoreGroup.get(i).getExtentInIndices();
            }
            
            // perform the scan by working your way back from the tail and by
            // constructing forward pointers as you go
            for(int i = coreCount - 2; i >= 0; i--)
            {
                List<IndexedSnpInterval> currCoreGroup = uberCores.get(i);
                int currCoreGroupSize = currCoreGroup.size();
                int[] currForwardPointers = new int[currCoreGroupSize];
                long[] currCumulativeExtents = new long[currCoreGroupSize];
                
                for(int j = 0; j < currCoreGroupSize; j++)
                {
                    IndexedSnpInterval currGroupInterval = currCoreGroup.get(j);
                    int prevCoreGroupSize = prevCoreGroup.size();
                    long maxCumulativeExtent = 0L;
                    for(int k = 0; k < prevCoreGroupSize; k++)
                    {
                        IndexedSnpInterval prevGroupInterval = prevCoreGroup.get(k);
                        long currCumulativeExtent =
                            cumulativeExtents[k] + currGroupInterval.getExtentInIndices();
                        if(currCumulativeExtent > maxCumulativeExtent &&
                           currGroupInterval.getEndIndex() >= prevGroupInterval.getStartIndex() - 1)
                        {
                            maxCumulativeExtent = currCumulativeExtent;
                            currCumulativeExtents[j] = currCumulativeExtent;
                            currForwardPointers[j] = k;
                        }
                    }
                    assert maxCumulativeExtent > 0L;
                }
                
                forwardPointers[i] = currForwardPointers;
                cumulativeExtents = currCumulativeExtents;
                prevCoreGroup = currCoreGroup;
            }
            
            // now we need to work our way forward through the pointers.
            // initialize the 1st pointer using the max cumulative extent
            // start with the max pointer
            int currPointer = 0;
            for(int i = 0; i < cumulativeExtents.length; i++)
            {
                if(cumulativeExtents[i] > cumulativeExtents[currPointer])
                {
                    currPointer = i;
                }
            }
            maxKIntervals.add(uberCores.get(0).get(currPointer));
            
            // the rest is easy. just hop through the forward pointers
            for(int i = 0; i < forwardPointers.length; i++)
            {
                currPointer = forwardPointers[i][currPointer];
                IndexedSnpInterval currMaxKInterval =
                    uberCores.get(i + 1).get(currPointer);
                maxKIntervals.add(currMaxKInterval);
            }
        }
        
        return maxKIntervals;
    }
    
    /**
     * Convert the indexed list into a physical position list. This will only
     * work if the intervals are ordered by thier starting positions
     * @param indexedSnpIntervals
     *        the index-based intervals
     * @param positionInputStream
     *        the position stream to do physical mapping
     * @return the intervals based on physical positions and
     * @throws IOException
     *         if we get an {@link IOException} while reading the stream
     */
    public List<BasePairInterval> toOrderedPhysicalIntervals(
            List<IndexedSnpInterval> indexedSnpIntervals,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        if(positionInputStream.getReadDirection() != StreamDirection.FORWARD)
        {
            throw new IllegalArgumentException(
                    "this function can only deal with forward reading " +
                    "position streams");
        }
        
        int chromosomeNumber = positionInputStream.getChromosomeNumber();
        
        int intervalCount = indexedSnpIntervals.size();
        List<BasePairInterval> physicalIntervals = new ArrayList<BasePairInterval>(
                intervalCount);
        if(intervalCount >= 1)
        {
            // initialize the position & index maps
            Map<Integer, Long> startIndexToStartPositionMap =
                new HashMap<Integer, Long>();
            Map<Integer, Integer> endIndexToStartIndexMap =
                new HashMap<Integer, Integer>();
            
            // for each indexed interval
            int positionStreamIndex = -1;
            long currPosition = -1;
            for(IndexedSnpInterval indexedInterval: indexedSnpIntervals)
            {
                int currStartIndex = indexedInterval.getStartIndex();
                assert startIndexToStartPositionMap.get(currStartIndex) == null;
                assert endIndexToStartIndexMap.get(indexedInterval.getEndIndex()) == null;
                assert indexedInterval.getExtentInIndices() >= 1;
                
                // move the position stream up to the same mark as the starting
                // index.
                while(positionStreamIndex < currStartIndex)
                {
                    // take the current position and increment the stream index
                    currPosition =
                        positionInputStream.getNextSnpPositionInBasePairs();
                    positionStreamIndex++;
                    
                    // does currPosition belong to one of our ending indices?
                    // if so we have completed a physical interval. lets
                    // remove it from the maps and add the interval to our list
                    Integer matchingStartIndex = endIndexToStartIndexMap.remove(
                            positionStreamIndex);
                    if(matchingStartIndex != null)
                    {
                        long matchingStartPosition =
                            startIndexToStartPositionMap.remove(matchingStartIndex);
                        physicalIntervals.add(new SimpleBasePairInterval(
                                chromosomeNumber,
                                matchingStartPosition,
                                1L + currPosition - matchingStartPosition));
                    }
                }
                assert positionStreamIndex == currStartIndex;
                
                if(indexedInterval.getExtentInIndices() == 1)
                {
                    physicalIntervals.add(new SimpleBasePairInterval(
                            chromosomeNumber,
                            currPosition,
                            1));
                }
                else
                {
                    startIndexToStartPositionMap.put(currStartIndex, currPosition);
                    endIndexToStartIndexMap.put(
                            indexedInterval.getEndIndex(),
                            currStartIndex);
                }
            }
            
            assert endIndexToStartIndexMap.size() == startIndexToStartPositionMap.size();
            
            // now handle the cleanup
            while(!endIndexToStartIndexMap.isEmpty())
            {
                currPosition =
                    positionInputStream.getNextSnpPositionInBasePairs();
                positionStreamIndex++;
                
                Integer matchingStartIndex = endIndexToStartIndexMap.remove(
                        positionStreamIndex);
                if(matchingStartIndex != null)
                {
                    long matchingStartPosition =
                        startIndexToStartPositionMap.remove(matchingStartIndex);
                    physicalIntervals.add(new SimpleBasePairInterval(
                            chromosomeNumber,
                            matchingStartPosition,
                            1L + currPosition - matchingStartPosition));
                }
            }
            assert startIndexToStartPositionMap.isEmpty();
        }
        assert physicalIntervals.size() == indexedSnpIntervals.size();
        
        // the intervals are ordered by ending index. we like them to
        // be ordered by starting index
        Collections.sort(physicalIntervals);
        return physicalIntervals;
    }
}
