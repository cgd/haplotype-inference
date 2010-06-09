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

package org.jax.haplotype.inference;

import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ArrayBlockingQueue;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;

/**
 * Similar in concept to {@link IntervalScanningHaplotypeEstimator} except the
 * underlying algorithms and interface are different.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SlidingWindowMultiHaplotypeEstimator implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 8295869216775611855L;

    private final ArrayBlockingQueue<SdpWithPosition> window;

    private final int windowSizeInSnps;

    private final boolean stepWindowBySnps;
    
    /**
     * A little class for holding on to both the SDP and the position of the
     * SDP
     */
    private class SdpWithPosition implements Serializable
    {
        /**
         * every {@link java.io.Serializable} is supposed to have one of these
         */
        private static final long serialVersionUID = -284658116323393931L;

        private final BitSet sdp;
        
        private final long position;

        /**
         * Constructor
         * @param sdp
         *          the SDP
         * @param position
         *          the position
         */
        public SdpWithPosition(BitSet sdp, long position)
        {
            this.sdp = sdp;
            this.position = position;
        }
        
        /**
         * Getter for the SDP
         * @return the sdp
         */
        public BitSet getSdp()
        {
            return this.sdp;
        }
        
        /**
         * Getter for the position
         * @return the position
         */
        public long getPosition()
        {
            return this.position;
        }
    }
    
    /**
     * Constructor
     * @param windowSizeInSnps
     *          the window size to use
     * @param stepWindowBySnps
     *          if true then the window will be stepped forward one snp at
     *          a time
     */
    public SlidingWindowMultiHaplotypeEstimator(
            int windowSizeInSnps,
            boolean stepWindowBySnps)
    {
        if(windowSizeInSnps <= 0)
        {
            throw new IllegalArgumentException(
                    "Window size must be greater positive not: " +
                    windowSizeInSnps);
        }
        
        this.windowSizeInSnps = windowSizeInSnps;
        this.stepWindowBySnps = stepWindowBySnps;
        
        this.window = new ArrayBlockingQueue<SdpWithPosition>(
                windowSizeInSnps,
                false);
    }
    
    /**
     * Getter for figuring out if we should step or slide the window.
     * if true then the window will be stepped forward one SNP at a time,
     * otherwise it is moved forward one whole window at a time
     * @return true if we're stepping
     */
    public boolean getStepWindowBySnps()
    {
        return this.stepWindowBySnps;
    }
    
    /**
     * Getter for the window size
     * @return
     *          the window size
     */
    public int getWindowSizeInSnps()
    {
        return this.windowSizeInSnps;
    }
    
    /**
     * Estimate the blocks. This function is NOT thread safe
     * @param sdpInputStream
     *          the SDP input stream
     * @param positionInputStream
     *          the position input stream
     * @return
     *          the blocks
     * @throws IOException
     *          if IO fails
     */
    public List<MultiPartitionedInterval> estimateMultiHaplotypeBlocks(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        ArrayList<MultiPartitionedInterval> hapBlocks =
            new ArrayList<MultiPartitionedInterval>();
        
        if(this.slideWindowForward(sdpInputStream, positionInputStream))
        {
//            long keepCount = 0;
//            long skipCount = 0;
            
            int chromosomeNumber = positionInputStream.getChromosomeNumber();
            int strainCount = sdpInputStream.getSdpStrainNames().length;
            MultiPartitionedInterval cumulativeBlock = null;
            do
            {
                MultiPartitionedInterval currBlock = this.estimateCurrentBlock(
                        chromosomeNumber,
                        strainCount);
                if(cumulativeBlock == null)
                {
                    cumulativeBlock = currBlock;
                }
                else
                {
                    if(Arrays.equals(
                            currBlock.getStrainGroups(),
                            cumulativeBlock.getStrainGroups()))
                    {
                        // since the strain groups match we should merge the two
                        // together into a single interval
                        long start = cumulativeBlock.getStartInBasePairs();
                        long end = currBlock.getEndInBasePairs();
                        long extent = 1 + end - start;
                        
                        cumulativeBlock = new MultiPartitionedInterval(
                                chromosomeNumber,
                                start,
                                extent,
                                cumulativeBlock.getStrainGroups());
                        
//                        skipCount++;
                    }
                    else
                    {
                        hapBlocks.add(cumulativeBlock);
                        cumulativeBlock = currBlock;
                        
//                        keepCount++;
                    }
                }
            } while(this.stepOrSlideWindowForward(sdpInputStream, positionInputStream));
            
            if(cumulativeBlock != null)
            {
                hapBlocks.add(cumulativeBlock);
                
//                keepCount++;
            }
            
//            System.out.println("Kept " + keepCount + " skipped " + skipCount);
        }
        
        hapBlocks.trimToSize();
        return hapBlocks;
    }
    
    /**
     * Estimate the block from the current window
     * @param strainCount
     *          the number of strains
     * @return
     *          the estimated partition
     */
    private MultiPartitionedInterval estimateCurrentBlock(
            int chromosomeNumber,
            int strainCount)
    {
        // start off by throwing every strain into the same group
        List<List<Integer>> strainGroups = new ArrayList<List<Integer>>();
        {
            List<Integer> everyStrain = new ArrayList<Integer>();
            for(int i = 0; i < strainCount; i++)
            {
                everyStrain.add(i);
            }
            strainGroups.add(everyStrain);
        }
        
        long startPos = -1;
        long endPos = -1;
        for(SdpWithPosition sdpWithPos: this.window)
        {
            // update the start and end positions
            if(startPos == -1)
            {
                startPos = sdpWithPos.getPosition();
            }
            endPos = sdpWithPos.getPosition();
            
            // for each SDP see if you can pull a group out that diverges
            // from the other strains in the group
            List<List<Integer>> divergentStrainGroups =
                new ArrayList<List<Integer>>();
            BitSet sdp = sdpWithPos.getSdp();
            for(List<Integer> strainGroup: strainGroups)
            {
                List<Integer> divergentStrainGroup = null;
                Iterator<Integer> strainGroupIter = strainGroup.iterator();
                
                // the first strain serves as a reference for the others
                boolean refBit = sdp.get(strainGroupIter.next());
                while(strainGroupIter.hasNext())
                {
                    int nextStrain = strainGroupIter.next();
                    if(sdp.get(nextStrain) != refBit)
                    {
                        // strain diverges from the reference, so remove it
                        // from this group and add it to the (possibly new)
                        // divergent group
                        strainGroupIter.remove();
                        if(divergentStrainGroup == null)
                        {
                            divergentStrainGroup = new ArrayList<Integer>();
                        }
                        
                        divergentStrainGroup.add(nextStrain);
                    }
                }
                
                if(divergentStrainGroup != null)
                {
                    divergentStrainGroups.add(divergentStrainGroup);
                }
            }
            
            strainGroups.addAll(divergentStrainGroups);
        }
        
        // this part just turns the grouping info into a more compact
        // representation (an array of shorts with length == # of strains and
        // value indicating group membership)
        short[] strainGroupings = new short[strainCount];
        short groupCount = (short)strainGroups.size();
        for(short groupNumber = 1; groupNumber < groupCount; groupNumber++)
        {
            List<Integer> currGroupMembers = strainGroups.get(groupNumber);
            for(int currGroupMember: currGroupMembers)
            {
                strainGroupings[currGroupMember] = groupNumber;
            }
        }
        
        return new MultiPartitionedInterval(
                chromosomeNumber,
                startPos,
                1 + endPos - startPos,
                strainGroupings);
    }

    /**
     * Either step or slide forward depending on {@link #getStepWindowBySnps()}
     * @param sdpInputStream
     *          the SDP input stream
     * @param positionInputStream
     *          the position input stream
     * @return
     *          true iff we haven't run out of data
     * @throws IOException
     *          if IO fails
     */
    private boolean stepOrSlideWindowForward(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        if(this.stepWindowBySnps)
        {
            return this.stepWindowForwardOneSnp(sdpInputStream, positionInputStream);
        }
        else
        {
            return this.slideWindowForward(sdpInputStream, positionInputStream);
        }
    }
    
    /**
     * Slide the window over {@link #getWindowSizeInSnps()} steps
     * @param sdpInputStream
     *          the SDP input stream
     * @param positionInputStream
     *          the position input stream
     * @return
     *          true if the fill succeeds
     * @throws IOException
     *          if IO fails
     */
    private boolean slideWindowForward(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        this.window.clear();
        for(int i = 0; i < this.windowSizeInSnps; i++)
        {
            if(!sdpInputStream.hasNextSdp() || !positionInputStream.hasNextSnpPosition())
            {
                return false;
            }
            
            this.window.add(new SdpWithPosition(
                    sdpInputStream.getNextSdp(),
                    positionInputStream.getNextSnpPositionInBasePairs()));
        }
        
        return true;
    }
    
    /**
     * Move the window forward one step. This function assumes that it is
     * moving an already full window (Eg from
     * {@link #slideWindowForward(SdpInputStream, SnpPositionInputStream)}) and
     * will leave a full window
     * @param sdpInputStream
     *          the SDP input stream
     * @param positionInputStream
     *          the position input stream
     * @return
     *          true unless we've reached the end of the streams
     * @throws IOException
     *          if IO fails
     */
    private boolean stepWindowForwardOneSnp(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException
    {
        this.window.remove();
        
        if(!sdpInputStream.hasNextSdp() || !positionInputStream.hasNextSnpPosition())
        {
            return false;
        }
        
        this.window.add(new SdpWithPosition(
                sdpInputStream.getNextSdp(),
                positionInputStream.getNextSnpPositionInBasePairs()));
        
        return true;
    }
}
