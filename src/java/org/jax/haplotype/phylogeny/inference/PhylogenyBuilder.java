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

package org.jax.haplotype.phylogeny.inference;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.SnpType;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdge;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.datastructure.SetUtilities;

/**
 * Class for building phylogenetic trees from
 * segments of the genome.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyBuilder
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            PhylogenyBuilder.class.getName());
    
    /**
     * Build phylogenetic trees
     * @param chromosomes
     *          the chromosomes to use
     * @param intervals
     *          the intervals to use
     * @return
     *          a list of phylogenetic trees. the list length will be equal
     *          to the number of intervals passed in
     * @throws NoValidPhylogenyException
     *          if we fail to generate any of the phylogeny trees
     */
    public List<PhylogenyInterval> buildPhylogeneticIntervals(
            List<StrainChromosome> chromosomes,
            List<BasePairInterval> intervals)
            throws NoValidPhylogenyException
    {
        if(!chromosomes.isEmpty() && !intervals.isEmpty())
        {
            SingleNucleotidePolymorphism[] anySnpList =
                chromosomes.get(0).getSingleNucleotidePolymorphisms();
            List<PhylogenyInterval> phylogenies = new ArrayList<PhylogenyInterval>();
            for(BasePairInterval interval: intervals)
            {
                // find the start index
                int startIndex = Arrays.binarySearch(
                        anySnpList,
                        new SingleNucleotidePolymorphism(
                                SnpType.A_SNP,
                                interval.getStartInBasePairs()),
                        SingleNucleotidePolymorphism.POSITION_ONLY_COMPARATOR);
                if(startIndex < 0)
                {
                    throw new IndexOutOfBoundsException(
                            "cant find position: " + interval.getStartInBasePairs());
                }
                
                // find the end index
                int endIndex = Arrays.binarySearch(
                        anySnpList,
                        new SingleNucleotidePolymorphism(
                                SnpType.A_SNP,
                                interval.getEndInBasePairs()),
                        SingleNucleotidePolymorphism.POSITION_ONLY_COMPARATOR);
                if(endIndex < 0)
                {
                    throw new IndexOutOfBoundsException(
                            "cant find position: " + interval.getEndInBasePairs());
                }
                
                phylogenies.add(this.buildPhylogeneticIntervals(
                        chromosomes,
                        new IndexedSnpInterval(
                                startIndex,
                                (endIndex - startIndex) + 1)));
            }
            return phylogenies;
        }
        else
        {
            return Collections.emptyList();
        }
    }

    /**
     * Build a phylogenetic tree
     * @param chromosomes
     *          the chromosomes to use for building the tree
     * @param interval
     *          the interval to generate a tree for
     * @return
     *          the tree
     * @throws NoValidPhylogenyException
     *          if we fail to generate the phylogeny tree
     */
    private PhylogenyInterval buildPhylogeneticIntervals(
            List<StrainChromosome> chromosomes,
            IndexedSnpInterval interval)
            throws NoValidPhylogenyException
    {
        // find all the unique SDPs in the interval
        int startIndex = interval.getStartIndex();
        int endIndex = interval.getEndIndex();
        Set<BitSet> sdpSet = new HashSet<BitSet>();
        List<BitSet> sdpList = new ArrayList<BitSet>();
        for(int snpIndex = startIndex; snpIndex <= endIndex; snpIndex++)
        {
            BitSet currSDP = this.convertSDPToBits(chromosomes, snpIndex);
            if(SetUtilities.isEmptySet(currSDP))
            {
                if(LOG.isLoggable(Level.FINE))
                {
                    LOG.fine("Ignoring SDP where all strains match");
                }
            }
            else if(sdpSet.add(currSDP))
            {
                sdpList.add(currSDP);
            }
        }
        
        // generate the tree from the SDPs
        PhylogenyTreeNode tree = this.buildPhylogeneticTreeFromSDPs(
                chromosomes,
                sdpList);
        
        assert new HashSet<String>(tree.getAllStrains()).size() ==
               tree.getAllStrains().size();
        assert tree.getAllStrains().size() == chromosomes.size();
        
        if(chromosomes.isEmpty())
        {
            return null;
        }
        else
        {
            BasePairInterval snpInterval =
                interval.toSnpInterval(chromosomes.get(0));
            return new PhylogenyInterval(tree, snpInterval);
        }
    }

    /**
     * Build the trees using the given SDP's using a version of
     * the algorithm presented by Ilan Gronau and Shlomo Moran in their
     * "Perfect Phylogeny" tutorial which was available at
     * http://webcourse.cs.technion.ac.il/236522/Spring2007/ho/WCFiles/tutorial10.ppt
     * @param chromosomes
     *          the chromosomes to build a tree from
     * @param sdpList
     *          the sdp list
     * @return
     *          the phylogeny tree
     * @throws NoValidPhylogenyException
     *          if there is no valid perfect phylogeny from the SDPs
     */
    private PhylogenyTreeNode buildPhylogeneticTreeFromSDPs(
            List<StrainChromosome> chromosomes,
            List<BitSet> sdpList) throws NoValidPhylogenyException
    {
        // sort the SDPs in an array from high to low
        List<BitSet> sortedSDPs = new ArrayList<BitSet>(sdpList);
        Collections.sort(
                sortedSDPs,
                Collections.reverseOrder(SetUtilities.BIT_SET_COMPARATOR));
        
        // recursively build 
        PhylogenyTreeNode rootPhylogeny = new PhylogenyTreeNode();
        while(!sortedSDPs.isEmpty())
        {
            BitSet childSDP = sortedSDPs.get(0);
            PhylogenyTreeNode childPhylogeny = this.buildPhylogeneticTreeFromSDPsRecursive(
                    chromosomes,
                    sortedSDPs,
                    0);
            rootPhylogeny.getChildEdges().add(new PhylogenyTreeEdge(
                    childSDP,
                    childPhylogeny,
                    1.0));
        }
        this.addStrainsToPhylogenyTreeNode(
                rootPhylogeny,
                null,
                chromosomes);
        return rootPhylogeny;
    }
    
    /**
     * Build the phylogenetic tree recursively. This recursive call has the
     * side effect of removing the SDP at SDP index in addition to any
     * child SDPs that fall under the edge of the SDP at the given index.
     * @see #buildPhylogeneticTreeFromSDPs(List, List)
     * @param chromosomes
     *          the chromosomes
     * @param sortedSDPs
     *          the sorted SDP values (this list is modified by this function)
     * @param sortedSDPIndex
     *          the SDP index to use on this recursive call
     * @return
     *          the tree node
     * @throws NoValidPhylogenyException
     *          if there is no valid perfect phylogeny from the SDPs
     */
    private PhylogenyTreeNode buildPhylogeneticTreeFromSDPsRecursive(
            List<StrainChromosome> chromosomes,
            List<BitSet> sortedSDPs,
            int sortedSDPIndex)
            throws NoValidPhylogenyException
    {
        // initialize current node's data
        PhylogenyTreeNode phylogeny = new PhylogenyTreeNode(
                new ArrayList<PhylogenyTreeEdge>(),
                new ArrayList<String>());
        BitSet sdpBits = sortedSDPs.get(sortedSDPIndex);
        
        // recursively add children
        {
            int i = sortedSDPIndex + 1;
            while(i < sortedSDPs.size())
            {
                BitSet nextSDP = sortedSDPs.get(i);
                if(sdpBits.intersects(nextSDP))
                {
                    // for perfect phylogeny the next SDP should be completely
                    // contained within the current SDP
                    if(!SetUtilities.isSubset(nextSDP, sdpBits))
                    {
                        throw new NoValidPhylogenyException(
                                "cannot build a perfect phylogeny with the " +
                                "following SDPs: " +
                                SetUtilities.bitSetToBinaryString(sdpBits) +
                                " and " +
                                SetUtilities.bitSetToBinaryString(nextSDP));
                    }
                    
                    // nextSDP is a child edge, recursively add it (note that
                    // we don't have do an i++ here because the children are
                    // recursively removed from the SDP list)
                    PhylogenyTreeNode childPhylogeny =
                        this.buildPhylogeneticTreeFromSDPsRecursive(
                                chromosomes,
                                sortedSDPs,
                                i);
                    PhylogenyTreeEdge childEdge = new PhylogenyTreeEdge(
                            nextSDP,
                            childPhylogeny,
                            1.0);
                    phylogeny.getChildEdges().add(childEdge);
                }
                else
                {
                    i++;
                }
            }
        }
        
        // add in the strains
        this.addStrainsToPhylogenyTreeNode(
                phylogeny,
                sdpBits,
                chromosomes);
        
        // we're done... remove us from the sorted SDP list so we don't get
        // added again
        sortedSDPs.remove(sortedSDPIndex);
        return phylogeny;
    }
    
    /**
     * Add the strains that live at this node. This function assumes that
     * the child edges have already been added with thier SDPs
     * @param phylogenyTreeNode
     *          the node to add strains to
     * @param nodeSDP
     *          the SDP of the edge that lead to this node. this should be
     *          null for the root node
     * @param chromosomes
     *          the list of chromosomes that match up with this bit ordering
     *          in the SDP
     */
    private void addStrainsToPhylogenyTreeNode(
            PhylogenyTreeNode phylogenyTreeNode,
            BitSet nodeSDP,
            List<StrainChromosome> chromosomes)
    {
        // we need to OR the children together for efficiency
        BitSet cumulativeChildrenSDP = new BitSet();
        for(PhylogenyTreeEdge childEdge: phylogenyTreeNode.getChildEdges())
        {
            cumulativeChildrenSDP.or(childEdge.getSdpBits());
        }
        
        // we need to figure out which strains (if any) live at this node. this
        // amounts to finding any strains that have a mutation (a '1' bit)
        // in this SDP but don't share mutations found in the SDP's of
        // child edges
        int strainCount = chromosomes.size();
        for(int i = 0; i < strainCount; i++)
        {
            if((nodeSDP == null || nodeSDP.get(i)) && !cumulativeChildrenSDP.get(i))
            {
                phylogenyTreeNode.getStrains().add(
                        chromosomes.get(i).getStrainName());
            }
        }
    }

    /**
     * Convert the strain distribution pattern (SDP) at the given snp
     * index into a bit set using a big integer to hold the bits where
     * 1 is the minority allele
     * @param chromosomes
     *          the chromosomes that we're calculating an SDP for
     * @param snpIndex
     *          the index to get the SDP at
     * @return
     *          the SDP
     */
    private BitSet convertSDPToBits(
            List<StrainChromosome> chromosomes,
            int snpIndex)
    {
        int chromosomeCount = chromosomes.size();
        
        if(chromosomeCount == 0)
        {
            throw new IllegalArgumentException(
                    "chromosome list cannot be empty");
        }
        else
        {
            // use the 1st strain as the reference
            BitSet sdpBits = new BitSet();
            SnpType referenceCall =
                chromosomes.get(0).getSingleNucleotidePolymorphisms()[snpIndex].getSnpType();
            
            int onesCount = 0;
            for(int i = 1; i < chromosomeCount; i++)
            {
                SnpType currCall =
                    chromosomes.get(i).getSingleNucleotidePolymorphisms()[snpIndex].getSnpType();
                if(currCall != referenceCall)
                {
                    sdpBits.set(i);
                    onesCount++;
                }
            }
            
            // make sure that 1 is the minority allele
            if(onesCount * 2 > chromosomeCount)
            {
                sdpBits.flip(0, chromosomeCount);
            }
            
            return sdpBits;
        }
    }
}
