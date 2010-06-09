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

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.ForwardStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.ForwardStrainChromosomeSnpPositionInputStream;
import org.jax.haplotype.io.GenotypeParser;
import org.jax.haplotype.io.ReverseStrainChromosomeSnpInputStream;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SimpleSdpInputStream;
import org.jax.haplotype.io.SnpInputStream;
import org.jax.util.datastructure.SequenceUtilities;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
@SuppressWarnings("all")
public class IntervalScannerTest
{
    @Test
    public void simpleScanTest() throws IOException
    {
        // Here's an example from Leonard's wiki page
        // SNP1    SNP2    SNP3    SNP4    SNP5    SNP6    SNP7    SNP8
        // A   0   0   0   0   0   0   0   0
        // B   1   0   0   0   1   0   1   1
        // C   1   1   0   0   0   0   1   0
        // D   0   0   1   0   0   1   1   1
        // E   0   0   0   1   0   1   1   1 
        String genoData =
            "snp.ID,ChrID,build.36.bp.Position,Source,A,B,C,D,E\n" +
            "snp1,1,1,test-src,T,G,G,T,T\n" +
            "snp2,1,2,test-src,T,T,G,T,T\n" +
            "snp3,1,3,test-src,G,G,G,T,G\n" +
            "snp4,1,4,test-src,G,G,G,G,T\n" +
            "snp5,1,5,test-src,G,T,G,G,G\n" +
            "snp6,1,6,test-src,G,G,G,T,T\n" +
            "snp7,1,7,test-src,G,T,T,T,T\n" +
            "snp8,1,8,test-src,G,A,G,A,A\n";
        
        GenotypeParser parser = new GenotypeParser();
        Set<StrainChromosome> chromosomes = parser.parseGenotypeFromStream(
                new ByteArrayInputStream(genoData.getBytes()));
        
        List<String> strainNamesList = parser.parseAvailableStrainsInOrderFromStream(
                new ByteArrayInputStream(genoData.getBytes()));
        
        String[] strainNames = strainNamesList.toArray(new String[0]);
        List<SnpInputStream> snpInputList = new ArrayList<SnpInputStream>();
        for(StrainChromosome strainChromosome: chromosomes)
        {
            ForwardStrainChromosomeSnpInputStream snpInputStream =
                new ForwardStrainChromosomeSnpInputStream(
                        chromosomes.iterator().next(),
                        strainChromosome);
            snpInputList.add(
                    snpInputStream);
        }
        SnpInputStream[] snpInputs = snpInputList.toArray(new SnpInputStream[0]);
        
        SdpInputStream sdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        IntervalScanner intervalScanner = new IntervalScanner();
        List<IndexedSnpInterval> intervals = intervalScanner.greedyScan(sdpInput);
        Assert.assertTrue(intervals.size() == 2);
    }
    
    @Test
    public void fullScanTestAllChromosomes() throws IOException
    {
        int chromosomes = 20;
        for(int i = 1; i <= chromosomes; i++)
        {
            System.out.println("=====================");
            System.out.println("FULL SCAN ITERATION " + i);
            this.fullScanTest("/chromosome_random_" + i + ".csv");
        }
    }
    
    public void fullScanTest(String genotypeResource) throws IOException
    {
        // initialization
        GenotypeParser parser = new GenotypeParser();
        Set<StrainChromosome> chromosomes = parser.parseGenotypeFromStream(
                IntervalScannerTest.class.getResourceAsStream(genotypeResource));
        List<String> strainNamesList = parser.parseAvailableStrainsInOrderFromStream(
                IntervalScannerTest.class.getResourceAsStream(genotypeResource));
        String[] strainNames = strainNamesList.toArray(new String[0]);
        ArrayList<SnpInputStream> snpInputList;
        SnpInputStream[] snpInputs;
        SdpInputStream sdpInput;
        IntervalScanner intervalScanner = new IntervalScanner();
        
        // Forward Scan
        StrainChromosome anyChromosome = chromosomes.iterator().next();
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
        sdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        System.out.println("Performing forward-scan");
        List<IndexedSnpInterval> forwardIntervals = intervalScanner.greedyScan(sdpInput);
        System.out.println("forward interval count: " + forwardIntervals.size());
//        Assert.assertTrue(forwardIntervals.size() == 11655);
        Assert.assertTrue(SequenceUtilities.isSorted(forwardIntervals));
        
        for(int i = 1; i < forwardIntervals.size(); i++)
        {
            Assert.assertTrue(
                    forwardIntervals.get(i - 1).getEndIndex() ==
                    forwardIntervals.get(i).getStartIndex() - 1);
        }
        
        // Reverse Scan
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
        sdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        System.out.println("Performing reverse-scan");
        List<IndexedSnpInterval> reverseIntervals = intervalScanner.greedyScan(sdpInput);
        System.out.println("reverse interval count: " + reverseIntervals.size());
        Assert.assertTrue(reverseIntervals.size() == forwardIntervals.size());
        Assert.assertTrue(SequenceUtilities.isSorted(reverseIntervals));
        
        for(int i = 1; i < reverseIntervals.size(); i++)
        {
            Assert.assertTrue(
                    reverseIntervals.get(i - 1).getEndIndex() ==
                    reverseIntervals.get(i).getStartIndex() - 1);
        }
        for(int i = 0; i < forwardIntervals.size(); i++)
        {
            IndexedSnpInterval forward = forwardIntervals.get(i);
            IndexedSnpInterval reverse = reverseIntervals.get(i);
            
            Assert.assertTrue(forward.getEndIndex() >= reverse.getEndIndex());
            Assert.assertTrue(reverse.getStartIndex() <= forward.getStartIndex());
        }
        
        // Uber Scan
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
        sdpInput = new SimpleSdpInputStream(
                strainNames,
                snpInputs);
        
        System.out.println("Performing uber-scan");
        List<IndexedSnpInterval> uberIntervals = intervalScanner.uberScan(sdpInput);
        System.out.println("uber interval count: " + uberIntervals.size());
//        Assert.assertTrue(uberIntervals.size() == 13393);
        Assert.assertTrue(SequenceUtilities.isSorted(uberIntervals));
        
        for(int i = 1; i < uberIntervals.size(); i++)
        {
            Assert.assertFalse(
                    uberIntervals.get(i - 1).contains(uberIntervals.get(i)));
            Assert.assertFalse(
                    uberIntervals.get(i).contains(uberIntervals.get(i - 1)));
        }
        
        // create cores
        List<IndexedSnpInterval> coreIntervals = intervalScanner.createCoreIntervals(
                forwardIntervals,
                reverseIntervals);
        Assert.assertTrue(coreIntervals.size() == forwardIntervals.size());
        Assert.assertTrue(SequenceUtilities.isSorted(coreIntervals));
        for(int i = 0; i < coreIntervals.size(); i++)
        {
            Assert.assertTrue(
                    forwardIntervals.get(i).contains(coreIntervals.get(i)));
            Assert.assertTrue(
                    reverseIntervals.get(i).contains(coreIntervals.get(i)));
            Assert.assertTrue(coreIntervals.get(i).getExtentInIndices() >= 1);
            Assert.assertTrue(
                    coreIntervals.get(i).getStartIndex() >= 0 &&
                    coreIntervals.get(i).getEndIndex() < sdpInput.getSdpCount());
            
            if(i >= 1)
            {
                Assert.assertTrue(
                        coreIntervals.get(i - 1).getEndIndex() <
                        coreIntervals.get(i).getStartIndex());
            }
        }
        
        // Subset uber intervals
        System.out.println("Subsetting uber intervals");
        List<List<IndexedSnpInterval>> uberCores = intervalScanner.createUberCores(
                uberIntervals,
                coreIntervals);
        Assert.assertTrue(uberCores.size() == forwardIntervals.size());
        for(int i = 0; i < uberCores.size(); i++)
        {
            List<IndexedSnpInterval> coreIntervalGroup = uberCores.get(i);
            
            Assert.assertTrue(
                    SequenceUtilities.isSorted(coreIntervalGroup));
            for(IndexedSnpInterval currCoreInterval: coreIntervalGroup)
            {
                Assert.assertTrue(currCoreInterval.contains(
                        coreIntervals.get(i)));
                if(i >= 1)
                {
                    Assert.assertFalse(currCoreInterval.intersects(
                            coreIntervals.get(i - 1)));
                }
                
                if(i + 1 < coreIntervals.size())
                {
                    Assert.assertFalse(currCoreInterval.intersects(
                            coreIntervals.get(i + 1)));
                }
            }
        }
        
        // Find maxk intervals
        System.out.println("Calculating MAX-K intervals");
        List<IndexedSnpInterval> maxKIntervals =
            intervalScanner.createMaxKIntervals(uberCores);
        Assert.assertTrue(maxKIntervals.size() == forwardIntervals.size());
        Assert.assertTrue(SequenceUtilities.isSorted(maxKIntervals));
        
        long indexedCumulativeExtent = 0;
        for(int i = 0; i < maxKIntervals.size(); i++)
        {
            indexedCumulativeExtent += maxKIntervals.get(i).getExtentInIndices();
            Assert.assertTrue(maxKIntervals.get(i).contains(coreIntervals.get(i)));
            
            if(i >= 1)
            {
                Assert.assertTrue(
                        maxKIntervals.get(i - 1).getEndIndex() >=
                        maxKIntervals.get(i).getStartIndex() - 1);
                Assert.assertFalse(
                        maxKIntervals.get(i - 1).contains(maxKIntervals.get(i)) ||
                        maxKIntervals.get(i).contains(maxKIntervals.get(i - 1)));
            }
        }
        
        long indexedTotalExtent = maxKIntervals.get(maxKIntervals.size() - 1).getEndIndex() + 1;
        Assert.assertTrue(indexedTotalExtent <= indexedCumulativeExtent);
        
        System.out.println("Indexed Max-k results:");
        System.out.println("Total Extent:       " + indexedTotalExtent);
        System.out.println("Cumulative Extent:  " + indexedCumulativeExtent);
        System.out.println("Cumulative Overlap: " + (indexedCumulativeExtent - indexedTotalExtent));
        
        // convert maxk to physical coordinates
        System.out.println("Converting max-k to physical coordinates");
        List<BasePairInterval> physicalMaxK = intervalScanner.toOrderedPhysicalIntervals(
                maxKIntervals,
                new ForwardStrainChromosomeSnpPositionInputStream(
                        anyChromosome));
        Assert.assertTrue(SequenceUtilities.isSorted(physicalMaxK));
        Assert.assertTrue(physicalMaxK.size() == forwardIntervals.size());
        
        long cumulativePhysicalExtent = 0L;
        for(int i = 0; i < physicalMaxK.size(); i++)
        {
            SingleNucleotidePolymorphism[] snps =
                anyChromosome.getSingleNucleotidePolymorphisms();
            int startIndex = maxKIntervals.get(i).getStartIndex();
            Assert.assertTrue(
                    snps[startIndex].getPositionInBasePairs() ==
                    physicalMaxK.get(i).getStartInBasePairs());
            int endIndex = maxKIntervals.get(i).getEndIndex();
            Assert.assertTrue(
                    snps[endIndex].getPositionInBasePairs() ==
                    physicalMaxK.get(i).getEndInBasePairs());
            cumulativePhysicalExtent += physicalMaxK.get(i).getExtentInBasePairs();
            if(i >= 1)
            {
                Assert.assertFalse(
                        physicalMaxK.get(i - 1).contains(physicalMaxK.get(i)));
                Assert.assertFalse(
                        physicalMaxK.get(i).contains(physicalMaxK.get(i - 1)));
            }
        }
        long totalPhysicalExtent =
            1 + physicalMaxK.get(physicalMaxK.size() - 1).getEndInBasePairs() -
            physicalMaxK.get(0).getStartInBasePairs();
        
        System.out.println("Physical Max-k Results");
        System.out.println("Total Extent:       " + totalPhysicalExtent);
        System.out.println("Cumulative Extent:  " + cumulativePhysicalExtent);
        System.out.println("Cumulative Overlap: " + (cumulativePhysicalExtent - totalPhysicalExtent));
        
        // Run the python version for comparison
//        System.out.println("running python code for validation");
//        CompgenPhylogenyInference phylogenyInference = new CompgenPhylogenyInference();
//        List<SimpleBasePairInterval> pyMaxKPhysicalIntervals = phylogenyInference.findMaxKCompatibleIntervals(
//                IntervalScannerTest.class.getResourceAsStream(genotypeResource));
//        Assert.assertTrue(forwardIntervals.size() == pyMaxKPhysicalIntervals.size());
//        System.out.println("Python Results:");
//        
//        long pyCumulativePhysicalExtent = 0L;
//        for(int i = 0; i < pyMaxKPhysicalIntervals.size(); i++)
//        {
//            pyCumulativePhysicalExtent += pyMaxKPhysicalIntervals.get(i).getExtentInBasePairs();
//        }
//        long pyTotalPhysicalExtent =
//            1 + pyMaxKPhysicalIntervals.get(pyMaxKPhysicalIntervals.size() - 1).getEndInBasePairs() -
//            pyMaxKPhysicalIntervals.get(0).getStartInBasePairs();
//        
//        System.out.println("Physical Python Max-k Results");
//        System.out.println("Total Extent:       " + pyTotalPhysicalExtent);
//        System.out.println("Cumulative Extent:  " + pyCumulativePhysicalExtent);
//        System.out.println("Cumulative Overlap: " + (pyCumulativePhysicalExtent - pyTotalPhysicalExtent));
//        
//        // convert python version to indices
//        List<IndexedSnpInterval> pyMaxKIntervals = this.toIndexedIntervals(
//                pyMaxKPhysicalIntervals,
//                anyChromosome.getSingleNucleotidePolymorphisms());
//        
//        long pyIndexedCumulativeExtent = 0;
//        for(int i = 0; i < pyMaxKIntervals.size(); i++)
//        {
//            pyIndexedCumulativeExtent += pyMaxKIntervals.get(i).getExtentInIndices();
//            Assert.assertTrue(pyMaxKIntervals.get(i).contains(coreIntervals.get(i)));
//            
//            if(i >= 1)
//            {
//                Assert.assertTrue(
//                        pyMaxKIntervals.get(i - 1).getEndIndex() >=
//                            pyMaxKIntervals.get(i).getStartIndex() - 1);
//                Assert.assertFalse(
//                        pyMaxKIntervals.get(i - 1).contains(pyMaxKIntervals.get(i)) ||
//                        pyMaxKIntervals.get(i).contains(pyMaxKIntervals.get(i - 1)));
//            }
//            
//            boolean matchInUberCores = false;
//            for(IndexedSnpInterval uberCoreInterval: uberCores.get(i))
//            {
//                if(uberCoreInterval.equals(pyMaxKIntervals.get(i)))
//                {
//                    matchInUberCores = true;
//                    break;
//                }
//            }
//            Assert.assertTrue(matchInUberCores);
//        }
//        
//        long pyIndexedTotalExtent = pyMaxKIntervals.get(pyMaxKIntervals.size() - 1).getEndIndex() + 1;
//        Assert.assertTrue(pyIndexedTotalExtent <= pyIndexedCumulativeExtent);
//        
//        System.out.println("Indexed Python Max-k results:");
//        System.out.println("Total Extent:       " + pyIndexedTotalExtent);
//        System.out.println("Cumulative Extent:  " + pyIndexedCumulativeExtent);
//        System.out.println("Cumulative Overlap: " + (pyIndexedCumulativeExtent - pyIndexedTotalExtent));
//        
//        Assert.assertTrue(indexedTotalExtent == pyIndexedTotalExtent);
//        Assert.assertTrue(indexedCumulativeExtent == pyIndexedCumulativeExtent);
    }

    private List<IndexedSnpInterval> toIndexedIntervals(
            List<BasePairInterval> physicalIntervals,
            SingleNucleotidePolymorphism[] snps)
    {
        List<IndexedSnpInterval> indexedIntervals = new ArrayList<IndexedSnpInterval>();
        for(BasePairInterval physicalInterval: physicalIntervals)
        {
            int startIndex = -1;
            for(int i = 0; i < snps.length; i++)
            {
                if(snps[i].getPositionInBasePairs() == physicalInterval.getStartInBasePairs())
                {
                    startIndex = i;
                }
                
                if(snps[i].getPositionInBasePairs() == physicalInterval.getEndInBasePairs())
                {
                    Assert.assertFalse(startIndex == -1);
                    indexedIntervals.add(new IndexedSnpInterval(
                            startIndex,
                            1 + i - startIndex));
                    break;
                }
            }
        }
        
        Assert.assertTrue(physicalIntervals.size() == indexedIntervals.size());
        for(int i = 0; i < indexedIntervals.size(); i++)
        {
            Assert.assertTrue(
                    snps[indexedIntervals.get(i).getStartIndex()].getPositionInBasePairs() ==
                    physicalIntervals.get(i).getStartInBasePairs());
            Assert.assertTrue(
                    snps[indexedIntervals.get(i).getEndIndex()].getPositionInBasePairs() ==
                    physicalIntervals.get(i).getEndInBasePairs());
        }
        
        return indexedIntervals;
    }
}
