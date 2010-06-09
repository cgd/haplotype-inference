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

package org.jax.haplotype.phylogeny.inference;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.geneticutil.data.SnpType;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.GenotypeParser;
import org.jax.haplotype.io.ForwardStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.ReverseStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SimpleSdpInputStream;
import org.jax.haplotype.io.SnpInputStream;
import org.jax.haplotype.phylogeny.data.NoValidPhylogenyException;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeEdge;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.datastructure.SequenceUtilities;
import org.junit.Assert;
import org.junit.Test;


/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
@SuppressWarnings("all")
public class PhylogenyScannerTest
{
    @Test
    public void phylogenyScanTestAllChromosomes() throws IOException, NoValidPhylogenyException
    {
        int chromosomes = 20;
        for(int i = 1; i <= chromosomes; i++)
        {
            System.out.println("=====================");
            System.out.println("PHYLOGENY SCAN ITERATION " + i);
            this.phylogenyScanTest("/chromosome_random_" + i + ".csv");
        }
    }

    private void phylogenyScanTest(String genotypeResource) throws IOException, NoValidPhylogenyException
    {
        // initialization
        GenotypeParser parser = new GenotypeParser();
        Set<StrainChromosome> chromosomeSet = parser.parseGenotypeFromStream(
                IntervalScannerTest.class.getResourceAsStream(genotypeResource));
        List<String> strainNamesList = parser.parseAvailableStrainsInOrderFromStream(
                IntervalScannerTest.class.getResourceAsStream(genotypeResource));
        String[] strainNames = strainNamesList.toArray(new String[0]);
        ArrayList<SnpInputStream> snpInputList;
        SnpInputStream[] snpInputs;
        SdpInputStream forwardSdpInput;
        SdpInputStream reverseSdpInput;
        SdpInputStream uberSdpInput;
        List<StrainChromosome> chromosomes = new ArrayList<StrainChromosome>();
        for(String strainName: strainNames)
        {
            StrainChromosome currOrderedChromosome = null;
            for(StrainChromosome currChromosome: chromosomeSet)
            {
                if(currChromosome.getStrainName().equals(strainName))
                {
                    currOrderedChromosome = currChromosome;
                }
            }
            Assert.assertTrue(currOrderedChromosome != null);
            chromosomes.add(currOrderedChromosome);
        }
        Assert.assertTrue(chromosomes.size() == strainNames.length);
        Assert.assertTrue(chromosomes.size() == chromosomeSet.size());
        
        StrainChromosome anyChromosome = chromosomes.iterator().next();
        IntervalScanner intervalScanner = new IntervalScanner();
        PhylogenyScanner phylogenyScanner = new PhylogenyScanner();
        
        // max-k Scan
        snpInputList = new ArrayList<SnpInputStream>();
        for(StrainChromosome strainChromosome: chromosomes)
        {
            ForwardStrainChromosomeSnpInputStream snpInputStream =
                new ForwardStrainChromosomeSnpInputStream(
                        anyChromosome,
                        strainChromosome);
            snpInputList.add(
                    snpInputStream);
        }
        snpInputs = snpInputList.toArray(new SnpInputStream[0]);
        forwardSdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        snpInputList = new ArrayList<SnpInputStream>();
        for(StrainChromosome strainChromosome: chromosomes)
        {
            ReverseStrainChromosomeSnpInputStream snpInputStream =
                new ReverseStrainChromosomeSnpInputStream(
                        anyChromosome,
                        strainChromosome);
            snpInputList.add(
                    snpInputStream);
        }
        snpInputs = snpInputList.toArray(new SnpInputStream[0]);
        reverseSdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        snpInputList = new ArrayList<SnpInputStream>();
        for(StrainChromosome strainChromosome: chromosomes)
        {
            ForwardStrainChromosomeSnpInputStream snpInputStream =
                new ForwardStrainChromosomeSnpInputStream(
                        anyChromosome,
                        strainChromosome);
            snpInputList.add(
                    snpInputStream);
        }
        snpInputs = snpInputList.toArray(new SnpInputStream[0]);
        uberSdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        System.out.println("Performing max k scan");
        List<IndexedSnpInterval> maxKIntervals = intervalScanner.maxKScan(
                forwardSdpInput,
                reverseSdpInput,
                uberSdpInput);
        System.out.println("done with max-k scan");
        
        System.out.println("forward interval count: " + maxKIntervals.size());
        Assert.assertTrue(SequenceUtilities.isSorted(maxKIntervals));
        
        snpInputList = new ArrayList<SnpInputStream>();
        for(StrainChromosome strainChromosome: chromosomes)
        {
            ForwardStrainChromosomeSnpInputStream snpInputStream =
                new ForwardStrainChromosomeSnpInputStream(
                        anyChromosome,
                        strainChromosome);
            snpInputList.add(
                    snpInputStream);
        }
        snpInputs = snpInputList.toArray(new SnpInputStream[0]);
        forwardSdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        System.out.println("inferring phylogenies");
        List<PhylogenyTreeNode> phylogenies = phylogenyScanner.inferPerfectPhylogenies(
                forwardSdpInput,
                maxKIntervals);
        System.out.println("done inferring phylogenies");
        Assert.assertTrue(phylogenies.size() == maxKIntervals.size());
        
//        CompgenPhylogenyInference compgenPhylogenyInference =
//            new CompgenPhylogenyInference();
//        List<SimpleBasePairInterval> compgenMaxkIntervals = compgenPhylogenyInference.findMaxKCompatibleIntervals(
//                PhylogenyScannerTest.class.getResourceAsStream(genotypeResource));
//        
//        assert compgenMaxkIntervals.size() == maxKIntervals.size();
//        
//        List<SimpleBasePairInterval> physicalMaxKIntervals = intervalScanner.toOrderedPhysicalIntervals(
//                maxKIntervals,
//                new ForwardStrainChromosomeSnpPositionInputStream(anyChromosome));
//        
//        List<PhylogenyInterval> phylogenyIntervals = compgenPhylogenyInference.buildPhylogeneticTrees(
//                PhylogenyScannerTest.class.getResourceAsStream(genotypeResource),
//                physicalMaxKIntervals);
//        
//        Assert.assertTrue(phylogenyIntervals.size() == phylogenies.size());
//        
//        for(int i = 0; i < phylogenies.size(); i++)
//        {
//            System.out.println(phylogenyIntervals.get(i).getPhylogeny().createNormalizedTree().toNewickFormat());
//            System.out.println(phylogenies.get(i).createNormalizedTree().toNewickFormat());
//            Assert.assertTrue(phylogenyIntervals.get(i).getPhylogeny().createNormalizedTree().equals(
//                    phylogenies.get(i).createNormalizedTree()));
//            if(i >= 1)
//            {
//                Assert.assertFalse(phylogenyIntervals.get(i).getPhylogeny().createNormalizedTree().equals(
//                        phylogenies.get(i - 1).createNormalizedTree()));
//            }
//        }
        
        List<StrainChromosome> chromosomeList = new ArrayList<StrainChromosome>();
        for(int i = 0; i < strainNames.length; i++)
        {
            for(StrainChromosome chromosome: chromosomes)
            {
                if(chromosome.getStrainName().equals(strainNames[i]))
                {
                    chromosomeList.add(chromosome);
                }
            }
        }
        
        Assert.assertTrue(chromosomeList.size() == strainNames.length);
        for(int i = 0; i < phylogenies.size(); i++)
        {
            PhylogenyTreeNode currPhylogeny = phylogenies.get(i);
            IndexedSnpInterval interval = maxKIntervals.get(i);
            
            Set<BitSet> phylogenySdps = this.calculateSdps(
                    currPhylogeny,
                    strainNames);
            Set<BitSet> intervalSdps = this.calculateSdps(
                    interval,
                    chromosomes,
                    strainNames);
            
            Assert.assertTrue(phylogenySdps.size() >= 1);
            Assert.assertTrue(phylogenySdps.size() <= interval.getExtentInIndices());
            Assert.assertTrue(intervalSdps.size() >= 1);
            Assert.assertTrue(intervalSdps.size() <= interval.getExtentInIndices());
//            System.out.println("phylogeny children: " + currPhylogeny.getChildEdges().size());
//            System.out.println("phylogeny SDPs: " + phylogenySdps.size());
//            for(BitSet sdp: phylogenySdps)
//            {
//                System.out.println(SetUtilities.bitSetToBinaryString(sdp));
//            }
//            System.out.println("interval SDPs:  " + intervalSdps.size());
//            for(BitSet sdp: intervalSdps)
//            {
//                System.out.println(SetUtilities.bitSetToBinaryString(sdp));
//            }
            // add the zero set to both
            phylogenySdps.add(new BitSet());
            intervalSdps.add(new BitSet());
            Assert.assertTrue(phylogenySdps.equals(intervalSdps));
//            Assert.assertTrue(phylogenySdps.containsAll(intervalSdps));
        }
    }
    
    private Set<BitSet> calculateSdps(
            IndexedSnpInterval interval,
            List<StrainChromosome> chromosomes,
            String[] strainNames) throws IOException
    {
//        ArrayList<SnpInputStream> snpInputList = new ArrayList<SnpInputStream>();
//        for(StrainChromosome strainChromosome: chromosomes)
//        {
//            ForwardStrainChromosomeSnpInputStream snpInputStream =
//                new ForwardStrainChromosomeSnpInputStream(
//                        chromosomes.get(0),
//                        strainChromosome);
//            snpInputList.add(
//                    snpInputStream);
//        }
//        SnpInputStream[] snpInputs = snpInputList.toArray(new SnpInputStream[0]);
//        SdpInputStream forwardSdpInput = new SimpleSdpInputStream(
//                strainNames,
//                snpInputs);
//        forwardSdpInput = new MinorityNormalizedSdpInputStream(forwardSdpInput);
//        
//        for(int i = 0; i < interval.getStartIndex(); i++)
//        {
//            forwardSdpInput.getNextSdp();
//        }
//        
//        Set<BitSet> sdps = new HashSet<BitSet>();
//        for(int i = 0; i < interval.getExtentInIndices(); i++)
//        {
//            sdps.add(forwardSdpInput.getNextSdp());
//        }
//        return sdps;
        Set<BitSet> sdps = new HashSet<BitSet>();
        for(int i = interval.getStartIndex(); i <= interval.getEndIndex(); i++)
        {
            BitSet sdp = PhylogenyScannerTest.convertSDPToBits(chromosomes, i);
            sdps.add(sdp);
        }
        
        return sdps;
    }
    
    private Set<BitSet> calculateSdps(
            PhylogenyTreeNode currPhylogeny,
            String[] allStrainNames)
    {
        Set<BitSet> sdpSet = new HashSet<BitSet>();
        
        
        for(PhylogenyTreeEdge childEdge: currPhylogeny.getChildEdges())
        {
            Set<BitSet> childSdpSet = this.calculateChildSdps(
                    childEdge.getNode(),
                    allStrainNames);
            
            // the bit set's from this child shouldn't intersect with
            // anything else that we've seen so far
            for(BitSet currChildSdp: childSdpSet)
            {
                for(BitSet otherChildSdp: sdpSet)
                {
                    Assert.assertFalse(
                            currChildSdp.intersects(otherChildSdp));
                }
            }
            
            sdpSet.addAll(childSdpSet);
        }
        
        BitSet cumulativeSdp = new BitSet();
        for(BitSet childSdp: sdpSet)
        {
            cumulativeSdp.or(childSdp);
        }
        
        List<String> allStrainNamesList = Arrays.asList(allStrainNames);
        for(String strainName: currPhylogeny.getStrains())
        {
            int index = allStrainNamesList.indexOf(strainName);
            
            Assert.assertTrue(index >= 0);
            Assert.assertFalse(cumulativeSdp.get(index));
        }
        
        return sdpSet;
    }

    private Set<BitSet> calculateChildSdps(
            PhylogenyTreeNode currPhylogeny,
            String[] allStrainNames)
    {
        Set<BitSet> sdpSet = new HashSet<BitSet>();
        
        for(PhylogenyTreeEdge childEdge: currPhylogeny.getChildEdges())
        {
            Set<BitSet> childSdpSet = this.calculateChildSdps(
                    childEdge.getNode(),
                    allStrainNames);
            
            // the bit set's from this child shouldn't intersect with
            // anything else that we've seen so far
            for(BitSet currChildSdp: childSdpSet)
            {
                for(BitSet otherChildSdp: sdpSet)
                {
                    Assert.assertFalse(
                            currChildSdp.intersects(otherChildSdp));
                }
            }
            
            sdpSet.addAll(childSdpSet);
        }
        
        BitSet currSdp = new BitSet();
        for(BitSet childSdp: sdpSet)
        {
            currSdp.or(childSdp);
        }
        
        List<String> allStrainNamesList = Arrays.asList(allStrainNames);
        for(String strainName: currPhylogeny.getStrains())
        {
            int index = allStrainNamesList.indexOf(strainName);
            
            Assert.assertTrue(index >= 0);
            Assert.assertFalse(currSdp.get(index));
            currSdp.set(index);
        }
        Assert.assertTrue(currSdp.cardinality() <= allStrainNames.length/2);
        
        sdpSet.add(currSdp);
        
        return sdpSet;
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
    public static BitSet convertSDPToBits(
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
