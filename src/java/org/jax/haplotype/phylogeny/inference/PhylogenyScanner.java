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
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.haplotype.io.MinorityNormalizedSdpInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdge;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.datastructure.SetUtilities;

/**
 * For building phylogenies using an SDP stream
 * @author Keith Sheppard
 */
public class PhylogenyScanner
{
    /**
     * Constructor
     */
    public PhylogenyScanner()
    {
    }
    
    /**
     * Infer perfect phylogenies using the given SDP inputs and the given
     * intervals. The intervals must be sorted for this function to work
     * @param sdpInputStream
     *          the SDP input stream
     * @param intervals
     *          the interval
     * @return
     *          the phylogenies with indices corresponding to the given
     *          intervals
     * @throws IOException
     *          if we get some kind of {@link IOException} reading the SDPs
     * @throws NoValidPhylogenyException
     *          if there isn't a valid phylogeny
     */
    public List<PhylogenyTreeNode> inferPerfectPhylogenies(
            SdpInputStream  sdpInputStream,
            List<IndexedSnpInterval> intervals) throws IOException, NoValidPhylogenyException
    {
        sdpInputStream = new MinorityNormalizedSdpInputStream(sdpInputStream);
        String[] sdpStrainNames = sdpInputStream.getSdpStrainNames();
        
        int intervalCount = intervals.size();
        
        Map<IndexedSnpInterval, PhylogenyTreeNode> intervalToPhylogenyMap =
            new HashMap<IndexedSnpInterval, PhylogenyTreeNode>();
        Map<IndexedSnpInterval, List<SdpInclusionHierarchy>> intervalToInclusionHierarchyMap =
            new HashMap<IndexedSnpInterval, List<SdpInclusionHierarchy>>();
        
        if(sdpInputStream.hasNextSdp())
        {
            int sdpIndex = 0;
            BitSet sdpBits = sdpInputStream.getNextSdp();
            for(IndexedSnpInterval interval: intervals)
            {
                int currStartIndex = interval.getStartIndex();
                assert sdpIndex <= currStartIndex;
                while(sdpIndex < currStartIndex)
                {
                    Iterator<Entry<IndexedSnpInterval, List<SdpInclusionHierarchy>>> hierarchyMapIter =
                        intervalToInclusionHierarchyMap.entrySet().iterator();
                    while(hierarchyMapIter.hasNext())
                    {
                        Entry<IndexedSnpInterval, List<SdpInclusionHierarchy>> currHierarchyMapEntry =
                            hierarchyMapIter.next();
                        IndexedSnpInterval currInterval =
                            currHierarchyMapEntry.getKey();
                        List<SdpInclusionHierarchy> inclusionHierarchies =
                            currHierarchyMapEntry.getValue();
                        
                        this.insertSdpInHierarchies(
                                inclusionHierarchies,
                                sdpBits);
                        if(currInterval.getEndIndex() <= sdpIndex)
                        {
                            assert currInterval.getEndIndex() == sdpIndex;
                            
                            // we're done with this interval
                            BitSet allBits = new BitSet(sdpStrainNames.length);
                            allBits.set(0, sdpStrainNames.length);
                            SdpInclusionHierarchy inclusionHierarchy =
                                new SdpInclusionHierarchy(
                                        allBits,
                                        inclusionHierarchies);
                            
                            PhylogenyTreeNode phylogeny = this.inclusionHierarchyToPhylogeny(
                                    inclusionHierarchy,
                                    sdpStrainNames);
                            assert !phylogeny.getChildEdges().isEmpty();
                            intervalToPhylogenyMap.put(
                                    currInterval,
                                    phylogeny);
                            hierarchyMapIter.remove();
                        }
                    }
                    
                    sdpBits = sdpInputStream.getNextSdp();
                    sdpIndex++;
                }
                
                intervalToInclusionHierarchyMap.put(
                        interval,
                        new ArrayList<SdpInclusionHierarchy>());
            }
            
            // clean up
            while(!intervalToInclusionHierarchyMap.isEmpty())
            {
                Iterator<Entry<IndexedSnpInterval, List<SdpInclusionHierarchy>>> hierarchyMapIter =
                    intervalToInclusionHierarchyMap.entrySet().iterator();
                while(hierarchyMapIter.hasNext())
                {
                    Entry<IndexedSnpInterval, List<SdpInclusionHierarchy>> currHierarchyMapEntry =
                        hierarchyMapIter.next();
                    IndexedSnpInterval currInterval =
                        currHierarchyMapEntry.getKey();
                    List<SdpInclusionHierarchy> inclusionHierarchies =
                        currHierarchyMapEntry.getValue();
                    
                    this.insertSdpInHierarchies(
                            inclusionHierarchies,
                            sdpBits);
                    if(currInterval.getEndIndex() <= sdpIndex)
                    {
                        assert currInterval.getEndIndex() == sdpIndex;
                        
                        // we're done with this interval
                        BitSet allBits = new BitSet(sdpStrainNames.length);
                        allBits.set(0, sdpStrainNames.length);
                        SdpInclusionHierarchy inclusionHierarchy =
                            new SdpInclusionHierarchy(
                                    allBits,
                                    inclusionHierarchies);
                        
                        PhylogenyTreeNode phylogeny = this.inclusionHierarchyToPhylogeny(
                                inclusionHierarchy,
                                sdpStrainNames);
                        assert !phylogeny.getChildEdges().isEmpty();
                        intervalToPhylogenyMap.put(
                                currInterval,
                                phylogeny);
                        hierarchyMapIter.remove();
                    }
                }
                
                if(sdpInputStream.hasNextSdp())
                {
                    sdpBits = sdpInputStream.getNextSdp();
                    sdpIndex++;
                }
                else
                {
                    assert intervalToInclusionHierarchyMap.isEmpty();
                }
            }
        }
        
        List<PhylogenyTreeNode> phylogenies =
            new ArrayList<PhylogenyTreeNode>(intervalCount);
        for(IndexedSnpInterval interval: intervals)
        {
            PhylogenyTreeNode currPhylogeny = intervalToPhylogenyMap.get(interval);
            assert currPhylogeny != null;
            phylogenies.add(currPhylogeny);
        }
        
        assert phylogenies.size() == intervalCount;
        return phylogenies;
    }
    
    /**
     * Convert an inclusion hierarchy into a phylogeny
     * @param inclusionHierarchy
     *          the inclusion hierarchy
     * @param sdpStrainNames
     *          the SDP strain names
     * @return
     *          the phylogeny node
     */
    private PhylogenyTreeNode inclusionHierarchyToPhylogeny(
            SdpInclusionHierarchy inclusionHierarchy,
            String[] sdpStrainNames)
    {
        // build the edge list
        ArrayList<PhylogenyTreeEdge> childEdges = new ArrayList<PhylogenyTreeEdge>(
                inclusionHierarchy.getChildren().size());
        BitSet combinedChildSdps = new BitSet();
        for(SdpInclusionHierarchy childInclusionHierarchy: inclusionHierarchy.getChildren())
        {
            combinedChildSdps.or(childInclusionHierarchy.getSdpBits());
            
            PhylogenyTreeNode childPhylogeny = this.inclusionHierarchyToPhylogeny(
                    childInclusionHierarchy,
                    sdpStrainNames);
            PhylogenyTreeEdge childEdge = new PhylogenyTreeEdge(
                    childInclusionHierarchy.getSdpBits(),
                    childPhylogeny,
                    1.0);
            childEdges.add(childEdge);
        }
        
        // build the strain name list
        ArrayList<String> strains = new ArrayList<String>();
        BitSet ourSdp = inclusionHierarchy.getSdpBits();
        for(int i = 0; i < sdpStrainNames.length; i++)
        {
            if(ourSdp.get(i) && !combinedChildSdps.get(i))
            {
                strains.add(sdpStrainNames[i]);
            }
        }
        
        // free up a little space
        childEdges.trimToSize();
        strains.trimToSize();
        
        return new PhylogenyTreeNode(
                childEdges,
                strains);
    }
    
    /**
     * Insert the SDP into a hierarchy
     * @param sdpHierarchyList
     *          the hierarchy list to modify
     * @param sdpBits
     *          the sdp bits
     * @throws NoValidPhylogenyException
     *          if the phylogeny isn't valid
     */
    private void insertSdpInHierarchies(
            List<SdpInclusionHierarchy> sdpHierarchyList,
            BitSet sdpBits) throws NoValidPhylogenyException
    {
        // nothing to do for empty sets
        if(!SetUtilities.isEmptySet(sdpBits))
        {
            this.insertSdpInHierarchiesRecursive(
                    sdpHierarchyList,
                    sdpBits);
        }
    }

    /**
     * Insert the SDP into a hierarchy recursively
     * @param sdpHierarchyList
     *          the hierarchy list to modify
     * @param sdpBits
     *          the sdp bits
     * @throws NoValidPhylogenyException
     *          if the phylogeny isn't valid
     */
    private void insertSdpInHierarchiesRecursive(
            List<SdpInclusionHierarchy> sdpHierarchyList,
            BitSet sdpBits) throws NoValidPhylogenyException
    {
        int sdpHierarchyListSize = sdpHierarchyList.size();
        for(int i = 0; i < sdpHierarchyListSize; i++)
        {
            SdpInclusionHierarchy currHierarchy = sdpHierarchyList.get(i);
            BitSet currSdpBits = currHierarchy.getSdpBits();
            if(currSdpBits.intersects(sdpBits))
            {
                if(currSdpBits.equals(sdpBits))
                {
                    // nothing to do. we've already accounted for the SDP
                    return;
                }
                else if(SetUtilities.isSubset(sdpBits, currSdpBits))
                {
                    // we're the current hierarchy's child
                    this.insertSdpInHierarchiesRecursive(
                            currHierarchy.getChildren(),
                            sdpBits);
                    return;
                }
                else if(SetUtilities.isSubset(currSdpBits, sdpBits))
                {
                    // the current hierarchy is our child
                    SdpInclusionHierarchy newSdpInclusionHierarchy =
                        new SdpInclusionHierarchy(sdpBits);
                    List<SdpInclusionHierarchy> newChildren =
                        newSdpInclusionHierarchy.getChildren();
                    newChildren.add(currHierarchy);
                    sdpHierarchyList.set(i, newSdpInclusionHierarchy);
                    
                    // see if we have any other children
                    for(int j = sdpHierarchyListSize - 1; j > i; j--)
                    {
                        currHierarchy = sdpHierarchyList.get(j);
                        currSdpBits = currHierarchy.getSdpBits();
                        
                        if(currSdpBits.intersects(sdpBits))
                        {
                            if(!SetUtilities.isSubset(currSdpBits, sdpBits))
                            {
                                throw new NoValidPhylogenyException(
                                        "can't create a perfect phylogeny. SDPs are " +
                                        "incompatible");
                            }
                            else
                            {
                                // move the child into the child list
                                // the remove is safe since we're counting
                                // down, not up
                                newChildren.add(currHierarchy);
                                sdpHierarchyList.remove(j);
                            }
                        }
                    }
                    
                    return;
                }
                else
                {
                    throw new NoValidPhylogenyException(
                            "can't create a perfect phylogeny. SDPs are " +
                            "incompatible");
                }
            }
        }
        
        SdpInclusionHierarchy newSdpInclusionHierarchy =
            new SdpInclusionHierarchy(sdpBits);
        sdpHierarchyList.add(newSdpInclusionHierarchy);
    }
}
