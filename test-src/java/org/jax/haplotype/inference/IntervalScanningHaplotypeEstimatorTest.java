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

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.EnumMap;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

import junit.framework.Assert;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.SnpType;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.CommaSeparatedChromosomeDataSource;
import org.jax.haplotype.inference.IntervalScanningHaplotypeEstimator.HaplotypeCandidate;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.util.SimpleTimer;
import org.jax.util.datastructure.SetUtilities;
import org.junit.Test;


/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntervalScanningHaplotypeEstimatorTest
{
    /**
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    @Test
    public void estimateHaplotypesTest() throws FileNotFoundException, IOException
    {
        Random random = new Random();
        
        IntervalScanningHaplotypeEstimator slidingWindowEstimator;
        for(int i = 1; i <= 20; i++)
        {
            ChromosomeDataSource chromoDataSource = new CommaSeparatedChromosomeDataSource(
                    IntervalScanningHaplotypeEstimatorTest.class.getResource("/chromosome_random_" + i + ".csv"),
                    i,
                    true);
            System.out.println("getting genotype");
            Set<StrainChromosome> genotype = chromoDataSource.getGenotypeData(
                    chromoDataSource.getAvailableStrains());
            System.out.println("done getting genotype");
            HashMap<String, StrainChromosome> genotypeMap = new HashMap<String, StrainChromosome>();
            for(StrainChromosome chromosome: genotype)
            {
                genotypeMap.put(chromosome.getStrainName(), chromosome);
            }
            
            String[] sortedStrains = chromoDataSource.getAvailableStrains().toArray(new String[0]);
            Arrays.sort(sortedStrains);
            System.out.println("getting SDPs");
            SdpInputStream sdpInput = chromoDataSource.getSdpInputStream(sortedStrains);
            System.out.println("getting positions");
            SnpPositionInputStream positionInput = chromoDataSource.getSnpPositionInputStream();
            System.out.println("done getting more stuff");

            int minNumSnps = random.nextInt(5) + 2;
            int minNumChromo = random.nextInt(2) + 2;
            slidingWindowEstimator =
                new IntervalScanningHaplotypeEstimator(
                        minNumSnps,
                        minNumChromo);
            
            SimpleTimer timer = new SimpleTimer();
            System.out.println("Estimating haplotypes");
            timer.reset();
            List<PartitionedInterval> estimatedHaplotypes = slidingWindowEstimator.estimateHaplotypeBlocks(
                    sdpInput,
                    positionInput);
            System.out.println("Time ellapsed: " + timer.getTimeEllapsedInSeconds() + "s");
            System.out.println("Sorting " + estimatedHaplotypes.size() + " haplotypes");
            Collections.sort(estimatedHaplotypes);
            
            System.out.println("Estimating haplotypes 2");
            timer.reset();
            List<PartitionedInterval> estimatedHaplotypes2 = this.oldEstimateHaplotypesFunction(
                    genotype,
                    minNumSnps,
                    minNumChromo);
            System.out.println("Time ellapsed: " + timer.getTimeEllapsedInSeconds() + "s");
            System.out.println("Sorting " + estimatedHaplotypes2.size() + " haplotypes");
            Collections.sort(estimatedHaplotypes2);
            Assert.assertEquals(estimatedHaplotypes.size(), estimatedHaplotypes2.size());
//            Assert.assertEquals(estimatedHaplotypes, estimatedHaplotypes2);
            for(int j = 0; j < estimatedHaplotypes.size(); j++)
            {
                if(!estimatedHaplotypes.get(j).equals(estimatedHaplotypes2.get(j)))
                {
                    System.out.println("Index " + j + " doesn't match");
                }
                Assert.assertEquals(
                        estimatedHaplotypes.get(j),
                        estimatedHaplotypes2.get(j));
            }
            
            System.out.println("==chromosome " + i + "==");
            System.out.println("number of snps:        " + genotype.iterator().next().getSingleNucleotidePolymorphisms().length);
            System.out.println("min hap snp count:     " + minNumSnps);
            System.out.println("min hap chromo group:  " + minNumChromo);
            System.out.println("num haps:              " + estimatedHaplotypes.size());
            
            System.out.print("Validating Haplotypes: ");
            Map<Long, Integer> positionToIndexMap = this.createPositionToIndexMap(
                    genotype);
            for(int j = 0; j < estimatedHaplotypes.size(); j += 1000)
            {
                int endingSnp = Math.min(estimatedHaplotypes.size(), j + 1000);
                this.checkHaplotypeRange(
                        sortedStrains,
                        genotypeMap,
                        positionToIndexMap,
                        estimatedHaplotypes,
                        j,
                        endingSnp,
                        minNumSnps,
                        minNumChromo);
                System.out.print("*");
            }
            System.out.println();
        }
    }
    
    /**
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    @Test
    public void lookForMissingHaplotypesTest() throws FileNotFoundException, IOException
    {
        Random random = new Random();
        
        IntervalScanningHaplotypeEstimator slidingWindowEstimator;
        String[] strainsToParse = new String[]{"129S1/SvImJ","129S4/SvJae", "129X1/SvJ","A/J","AKR/J","BALB/cByJ","BUB/BnJ","C3H/HeJ","C57BL/6J","C57BLKS/J","C57BR/cdJ","C57L/J","C58/J","CAST/EiJ","CBA/J","CE/J","CZECHII/EiJ","DBA/1J"};
        Arrays.sort(strainsToParse);
        
        for(int iteration = 0; iteration < 1; iteration++)
        {
            System.out.println("Look for missing iteration: " + iteration);
            for(int i = 1; i <= 20; i++)
//            for(int i = 5; i <= 10; i++)
            {
                ChromosomeDataSource chromoDataSource = new CommaSeparatedChromosomeDataSource(
                        IntervalScanningHaplotypeEstimatorTest.class.getResource("/chromosome_random_" + i + ".csv"),
                        i,
                        true);
                Assert.assertTrue(
                        chromoDataSource.getAvailableStrains().containsAll(
                        Arrays.asList(strainsToParse)));
                int minNumSnps = 1 + random.nextInt(3);
                int minNumChromo = 2 + random.nextInt(3);
                slidingWindowEstimator =
                    new IntervalScanningHaplotypeEstimator(
                            minNumSnps,
                            minNumChromo);
                System.out.println("Estimating haplotypes");
                List<PartitionedInterval> estimatedHaplotypes = slidingWindowEstimator.estimateHaplotypeBlocks(
                        chromoDataSource.getSdpInputStream(strainsToParse),
                        chromoDataSource.getSnpPositionInputStream());
                System.out.println("Sorting haplotypes: " + estimatedHaplotypes.size());
                Collections.sort(estimatedHaplotypes);
                System.out.println("Estimating haplotypes 2");
                List<PartitionedInterval> estimatedHaplotypes2 = this.oldEstimateHaplotypesFunction(
                        chromoDataSource.getGenotypeData(new HashSet<String>(Arrays.asList(strainsToParse))),
                        minNumSnps,
                        minNumChromo);
                System.out.println("Sorting haplotypes 2: " + estimatedHaplotypes2.size());
                Collections.sort(estimatedHaplotypes2);
                
//                Assert.assertEquals(estimatedHaplotypes, estimatedHaplotypes2);
                
                System.out.println("==chromosome " + i + "==");
                System.out.println("number of snps:        " + chromoDataSource.getSnpPositionInputStream().getSnpCount());
                System.out.println("min hap snp count:     " + minNumSnps);
                System.out.println("min hap chromo group:  " + minNumChromo);
                System.out.println("num haps:              " + estimatedHaplotypes.size());
                
                Set<StrainChromosome> genotype = chromoDataSource.getGenotypeData(
                        new HashSet<String>(Arrays.asList(strainsToParse)));
                StrainChromosome[] genotypeArray = genotype.toArray(new StrainChromosome[0]);
                Arrays.sort(genotypeArray);
                
                Assert.assertTrue(genotypeArray.length == strainsToParse.length);
                for(int j = 0; j < genotypeArray.length; j++)
                {
                    Assert.assertEquals(genotypeArray[j].getStrainName(), strainsToParse[j]);
                }
                
                this.lookForMissingHaplotypes(
                        genotypeArray,
                        estimatedHaplotypes2,
                        minNumSnps,
                        minNumChromo,
                        random);
                this.lookForMissingHaplotypes(
                        genotypeArray,
                        estimatedHaplotypes,
                        minNumSnps,
                        minNumChromo,
                        random);
                
                System.out.println();
            }
        }
    }
    
    private void lookForMissingHaplotypes(
            StrainChromosome[] genotypeArray,
            List<PartitionedInterval> estimatedHaplotypes,
            int minimumExtent,
            int minimumGroupSize,
            Random random)
    {
        System.out.println("Looking for missing haplotypes");
//        for(int i = 0; i < 10; i++)
        {
            int randomIndex1 = random.nextInt(genotypeArray.length);
            
            int randomIndex2 = random.nextInt(genotypeArray.length - 1);
            if(randomIndex2 >= randomIndex1)
            {
                randomIndex2++;
            }
            
            System.out.println("random strain index 1: " + randomIndex1);
            System.out.println("random strain 1: " + genotypeArray[randomIndex1].getStrainName());
            System.out.println("random strain index 2: " + randomIndex2);
            System.out.println("random strain 2: " + genotypeArray[randomIndex2].getStrainName());
            
            StrainChromosome chromo1 = genotypeArray[randomIndex1];
            StrainChromosome chromo2 = genotypeArray[randomIndex2];
            SingleNucleotidePolymorphism[] snps1 = chromo1.getSingleNucleotidePolymorphisms();
            SingleNucleotidePolymorphism[] snps2 = chromo2.getSingleNucleotidePolymorphisms();
            
            assert snps1.length == snps2.length;
            
            for(int currSnp = 0; currSnp < snps1.length; currSnp++)
            {
                int agreementStart = currSnp;
                while(currSnp < snps1.length && snps1[currSnp].getSnpType() == snps2[currSnp].getSnpType())
                {
                    currSnp++;
                }
                
                if(currSnp - agreementStart >= minimumExtent)
                {
                    BitSet strainsInHaplotype = null;
                    for(int agreementSnpCursor = agreementStart;
                       agreementSnpCursor < currSnp;
                       agreementSnpCursor++)
                    {
                        BitSet[] groups = IntervalScanningHaplotypeEstimatorTest.groupChromosomesBySnpTypeAt(
                                genotypeArray,
                                agreementSnpCursor);
                        
                        BitSet matchingGroup = null;
                        for(BitSet group: groups)
                        {
//                            List<String> groupAsList = Arrays.asList(group);
//                            if(groupAsList.contains(chromo1.getStrainName()))
//                            {
//                                assert groupAsList.contains(chromo2.getStrainName());
//                                assert matchingGroupList == null;
//                                matchingGroupList = groupAsList;
//                            }
                            if(group.get(randomIndex1))
                            {
                                assert group.get(randomIndex2);
                                assert matchingGroup == null;
                                matchingGroup = group;
                            }
                        }
                        
                        if(strainsInHaplotype == null)
                        {
                            strainsInHaplotype = matchingGroup;
                        }
                        else
                        {
                            strainsInHaplotype.and(matchingGroup);
                        }
                    }
                    
                    assert strainsInHaplotype.cardinality() >= 2;
                    
                    if(strainsInHaplotype.cardinality() >= minimumGroupSize)
                    {
//                        System.out.print('+');
                        long startInBases = snps1[agreementStart].getPositionInBasePairs();
                        long endInBases = snps1[currSnp - 1].getPositionInBasePairs();
                        
                        boolean foundMatch = false;
                        for(PartitionedInterval haplotypeBlock: estimatedHaplotypes)
                        {
                            if(haplotypeBlock.getStartInBasePairs() == startInBases)
                            {
                                if(haplotypeBlock.getEndInBasePairs() == endInBases)
                                {
                                    if(haplotypeBlock.getStrainBitSet().get(randomIndex1) &&
                                       haplotypeBlock.getStrainBitSet().get(randomIndex2))
                                    {
    //                                    System.out.println("match for: " + startInBases + " " + endInBases);
                                        foundMatch = true;
                                        break;
                                    }
                                }
                            }
                        }
                        
                        if(!foundMatch)
                        {
                            String errorMessage =
                                "Missed Haplotype: " + chromo1.getStrainName() +
                                " " + chromo2.getStrainName() + " start bases:" +
                                startInBases + " end bases: " + endInBases +
                                " start index: " + agreementStart + " size snps: " +
                                (currSnp - agreementStart);
                            System.out.println(errorMessage);
                            
//                            System.out.println("inferred haplotypes with matching range:");
//                            for(PartitionedInterval haplotypeBlock: estimatedHaplotypes)
//                            {
//                                if(haplotypeBlock.getStartInBasePairs() == startInBases)
//                                {
//                                    if(haplotypeBlock.getEndingSnpInBasePairs() == endInBases)
//                                    {
//                                        System.out.println("Haplotype Strain Group:");
//                                        for(String strainName: haplotypeBlock.getStrainBitSet())
//                                        {
//                                            System.out.println(strainName);
//                                        }
//                                    }
//                                }
//                            }
                            
                            System.out.println("all strains that should be in haplotype:");
                            System.out.println(
                                    SetUtilities.bitSetToBinaryString(
                                            strainsInHaplotype));
    
                            throw new RuntimeException(
                                    errorMessage);
                        }
                    }
//                    else
//                    {
//                        System.out.print('-');
//                    }
                }
            }
//            System.out.print("*");
        }
//        System.out.println();
        
        System.out.println("Didn't find missing haplotypes");
    }

    /**
     * 
     * @throws FileNotFoundException
     * @throws IOException
     */
    @Test
    public void equivalenceClassTest() throws FileNotFoundException, IOException
    {
        int minNumSnps = 3;
        int minNumChromo = 5;
        IntervalScanningHaplotypeEstimator slidingWindowEstimator;
        List<PartitionedInterval> allHaplotypeBlocks =
            new ArrayList<PartitionedInterval>();
        int cumulativeChromosomeEquivClasses = 0;
        for(int i = 1; i <= 20; i++)
        {
            String resourceName =
                "/chromosome_random_" + i + ".csv";
            ChromosomeDataSource chromoData = new CommaSeparatedChromosomeDataSource(
                    IntervalScanningHaplotypeEstimatorTest.class.getResource(resourceName),
                    i,
                    true);
            
//            Set<StrainChromosome> genotype = parser.parseGenotypeFromStream(
////                    new FileInputStream(file),
//                    IntervalScanningHaplotypeEstimatorTest.class.getResourceAsStream(resourceName),
//                    true);
//            HashMap<String, StrainChromosome> genotypeMap = new HashMap<String, StrainChromosome>();
//            for(StrainChromosome chromosome: genotype)
//            {
//                genotypeMap.put(chromosome.getStrainName(), chromosome);
//            }
            System.out.println("estimating haps");
            slidingWindowEstimator =
                new IntervalScanningHaplotypeEstimator(
                        minNumSnps,
                        minNumChromo);
            String[] strains = chromoData.getAvailableStrains().toArray(new String[0]);
            Arrays.sort(strains);
            List<PartitionedInterval> estimatedHaplotypes = slidingWindowEstimator.estimateHaplotypeBlocks(
                    chromoData.getSdpInputStream(strains),
                    chromoData.getSnpPositionInputStream());
            Collections.sort(estimatedHaplotypes);
            
            System.out.println(
                    "Number of haplotype blocks for chromosome " + i + ": " +
                    estimatedHaplotypes.size());
            List<PartitionedIntervalSet> currEquivalenceClasses =
                HaplotypeEquivalenceClassCreator.createEquivalenceClassesFromBlocks(
                        estimatedHaplotypes);
            System.out.println(
                    "Number of unique strain groups for chromosome " + i +
                    ": " + currEquivalenceClasses.size());
            cumulativeChromosomeEquivClasses += currEquivalenceClasses.size();
            
//            for(int j = 0; j < estimatedHaplotypes.size(); j += 1000)
//            {
//                int endingSnp = Math.min(estimatedHaplotypes.size(), j + 1000);
//                this.checkHaplotypeRange(
//                        genotypeMap,
//                        estimatedHaplotypes,
//                        j,
//                        endingSnp,
//                        minNumSnps,
//                        minNumChromo);
//            }
            allHaplotypeBlocks.addAll(estimatedHaplotypes);
        }
        
        List<PartitionedIntervalSet> allEquivalenceClasses =
            HaplotypeEquivalenceClassCreator.createEquivalenceClassesFromBlocks(
                allHaplotypeBlocks);
        System.out.println(
                "total number of equivalence classes: " + allEquivalenceClasses.size());
        System.out.println(
                "cumulative chromosome only equivalence classes: " +
                cumulativeChromosomeEquivClasses);
    }
    
    private Map<Long, Integer> createPositionToIndexMap(
            Set<StrainChromosome> genotype)
    {
        StrainChromosome chromo = genotype.iterator().next();
        SingleNucleotidePolymorphism[] snps = chromo.getSingleNucleotidePolymorphisms();
        HashMap<Long, Integer> positionToIndexMap = new HashMap<Long, Integer>(
                snps.length);
        for(int i = 0; i < snps.length; i++)
        {
            SingleNucleotidePolymorphism snp = snps[i];
            positionToIndexMap.put(snp.getPositionInBasePairs(), i);
        }
        
        return positionToIndexMap;
    }

    private void checkHaplotypeRange(
            String[] sortedStrainNames,
            Map<String, StrainChromosome> chromosomes,
            Map<Long, Integer> positionToIndexMap,
            List<PartitionedInterval> haplotypeList,
            int startingSnp,
            int endingSnpExclusive,
            int windowSize,
            int minNumChromosomes)
    {
        for(int i = startingSnp; i < endingSnpExclusive; i++)
        {
            PartitionedInterval haplotype1 = haplotypeList.get(i);
            this.validateHaplotype(sortedStrainNames, chromosomes, positionToIndexMap, haplotype1, windowSize, minNumChromosomes);
            
            for(int j = startingSnp; j < endingSnpExclusive; j++)
            {
                PartitionedInterval haplotype2 = haplotypeList.get(j);
                if(i != j)
                {
                    this.checkForHaplotypeConflicts(haplotype1, haplotype2);
                }
            }
        }
    }
    
    /**
     * @param haplotype
     * @param haplotype2
     */
    private void checkForHaplotypeConflicts(
            PartitionedInterval haplotype,
            PartitionedInterval haplotype2)
    {
        IntervalScanningHaplotypeEstimatorTest.checkHaplotype1OverlapsHaplotype2(haplotype, haplotype2);
    }
    
    private static void checkHaplotype1OverlapsHaplotype2(
            PartitionedInterval haplotype1,
            PartitionedInterval haplotype2)
    {
        if(haplotype1.getStartInBasePairs() <= haplotype2.getStartInBasePairs() &&
           haplotype1.getEndInBasePairs() >= haplotype2.getEndInBasePairs())
        {
            BitSet intersection = (BitSet)haplotype1.getStrainBitSet().clone();
            intersection.and(haplotype2.getStrainBitSet());
            
            if(intersection.cardinality() == haplotype2.getStrainBitSet().cardinality())
            {
                throw new IllegalArgumentException(
                       "haplotype 1 overlaps haplotype2");
            }
        }
    }
    
    /**
     * @param currHap
     */
    private void validateHaplotype(
            String[] sortedStrains,
            Map<String, StrainChromosome> chromosomes,
            Map<Long, Integer> positionToIndexMap,
            PartitionedInterval currHap,
            int windowSize,
            int minNumChromo)
    {
        if(currHap.getStrainBitSet().cardinality() < minNumChromo)
        {
            throw new IllegalArgumentException(
                    "too few chromosomes in haplotype: " +
                    currHap.getStrainBitSet().cardinality());
        }
        
        if(currHap.getExtentInBasePairs() < windowSize)
        {
            throw new IllegalArgumentException(
                    "too few snps " + currHap.getExtentInBasePairs());
        }
        
        int startingSnpIndex = positionToIndexMap.get(currHap.getStartInBasePairs());
        int endingSnpIndex = positionToIndexMap.get(currHap.getEndInBasePairs());
        int[] haplotypeStrainIndeces = SetUtilities.getSetBitIndices(currHap.getStrainBitSet());
        for(int currSnpIndex = startingSnpIndex; currSnpIndex <= endingSnpIndex; currSnpIndex++)
        {
            for(int j = 1; j < haplotypeStrainIndeces.length; j++)
            {
                StrainChromosome currChromosome = chromosomes.get(
                    sortedStrains[haplotypeStrainIndeces[j]]);
                StrainChromosome prevChromosome = chromosomes.get(
                        sortedStrains[haplotypeStrainIndeces[j - 1]]);
                if(currChromosome.getSingleNucleotidePolymorphisms()[currSnpIndex].getSnpType() !=
                   prevChromosome.getSingleNucleotidePolymorphisms()[currSnpIndex].getSnpType())
                {
                    throw new IllegalArgumentException(
                            "snps dont match at index " + currSnpIndex +
                            " and chromosome " + j);
                }
            }
        }
    }
    
    private List<PartitionedInterval> oldEstimateHaplotypesFunction(
            Set<StrainChromosome> chromosomes,
            int windowSizeInSnps,
            int minimumNumberOfChromosomes)
    {
        List<PartitionedInterval> haplotypeList =
            new ArrayList<PartitionedInterval>();
        StrainChromosome[] sortedChromosomes =
            chromosomes.toArray(new StrainChromosome[chromosomes.size()]);
        Arrays.sort(sortedChromosomes);
        
        IntervalScanningHaplotypeEstimatorTest.validateChromosomes(sortedChromosomes);
        
        if(sortedChromosomes.length >= 2)
        {
            int snpCount =
                sortedChromosomes[0].getSingleNucleotidePolymorphisms().length;
            
            List<HaplotypeCandidate> haplotypeCandidateList =
                new ArrayList<HaplotypeCandidate>();
            List<HaplotypeCandidate> newHaplotypeCandidateList =
                new ArrayList<HaplotypeCandidate>();
            List<HaplotypeCandidate> potentialNewCandidateList =
                new ArrayList<HaplotypeCandidate>();
            
            // iterate through the SNP positions
            for(int currSnpIndex = 0; currSnpIndex < snpCount; currSnpIndex++)
            {
                BitSet[] snpGroups =
                    IntervalScanningHaplotypeEstimatorTest.groupChromosomesBySnpTypeAt(
                            sortedChromosomes,
                            currSnpIndex);
                assert snpGroups.length >= 1;
                
                // reinitialize lists
                newHaplotypeCandidateList.clear();
                potentialNewCandidateList.clear();
                
                // add snp groups as candidates
                for(BitSet currSnpGroup: snpGroups)
                {
                    if(currSnpGroup.cardinality() >= minimumNumberOfChromosomes)
                    {
                        potentialNewCandidateList.add(new HaplotypeCandidate(
                                -1, // don't care
                                currSnpIndex,
                                currSnpGroup));
                    }
                }
                
                for(HaplotypeCandidate currHaplotypeCandidate: haplotypeCandidateList)
                {
                    boolean terminateCandidateHaplotype = false;
                    
                    BitSet currHaplotypeCandidateStrains =
                        currHaplotypeCandidate.getStrainBitSet();
                    int currHaplotypeCandidateBitCount =
                        currHaplotypeCandidateStrains.cardinality();
                    SNP_GROUP_LOOP:
                    for(BitSet currSnpGroup: snpGroups)
                    {
                        // create an intersection of the current snp group and the
                        // current candidate
                        BitSet chromosomeIntersection = (BitSet)currSnpGroup.clone();
                        chromosomeIntersection.and(currHaplotypeCandidateStrains);
                        
                        int intersectionBitCount = chromosomeIntersection.cardinality();
                        assert intersectionBitCount <= currSnpGroup.cardinality();
                        assert intersectionBitCount <= currHaplotypeCandidateBitCount;
                        
                        if(intersectionBitCount > 0)
                        {
                            // see if all of the haplotype chromosomes are
                            // contained in the current SNP group
                            if(currHaplotypeCandidateBitCount ==
                               intersectionBitCount)
                            {
                                assert !terminateCandidateHaplotype;
                                
                                // Pass the haplotype through to the new
                                // candidate list unchanged
                                newHaplotypeCandidateList.add(
                                        currHaplotypeCandidate);
                                
                                // no other SNP groups will overlap with this
                                // haplotype candidate
                                break SNP_GROUP_LOOP;
                            }
                            else
                            {
                                // it's a partial intersection.
                                // 1st terminate the old potential haplotype
                                // since it's been broken up
                                terminateCandidateHaplotype = true;
                                
                                // 2nd add the new haplotype to the candidate
                                // list if it is big enough
                                if(intersectionBitCount >= minimumNumberOfChromosomes)
                                {
                                    // this new potential candidate should go back
                                    // as far as the current candidate does
                                    potentialNewCandidateList.add(new HaplotypeCandidate(
                                            -1, // don't care
                                            currHaplotypeCandidate.getStartingSnpIndex(),
                                            chromosomeIntersection));
                                }
                            }
                        }
                    }
                    
                    // see if the haplotype was terminated
                    if(terminateCandidateHaplotype)
                    {
                        int indexDifference =
                            currSnpIndex - currHaplotypeCandidate.getStartingSnpIndex();
                        if(indexDifference >= windowSizeInSnps)
                        {
                            // the window size is big enough to put
                            // into the permanent haplotype list
                            haplotypeList.add(IntervalScanningHaplotypeEstimatorTest.createHaplotypeBlockFromCandidate(
                                    currHaplotypeCandidate,
                                    currSnpIndex,
                                    sortedChromosomes));
                        }
                    }
                }
                
                // add potential candidates only if they pass through the filter
                // filter against other potential candidates 1st
                {
                    Iterator<HaplotypeCandidate> potentialCandidateIter =
                        potentialNewCandidateList.iterator();
                    while(potentialCandidateIter.hasNext())
                    {
                        HaplotypeCandidate currPotentialCandidate =
                            potentialCandidateIter.next();
                        
                        for(HaplotypeCandidate otherPotentialCandidate: potentialNewCandidateList)
                        {
                            if(currPotentialCandidate != otherPotentialCandidate)
                            {
                                if(this.doesCandidate1ContainCandidate2(
                                        otherPotentialCandidate,
                                        currPotentialCandidate))
                                {
                                    assert !currPotentialCandidate.equals(otherPotentialCandidate);
                                    
                                    // remove this candidate since its fully contained
                                    potentialCandidateIter.remove();
                                    break;
                                }
                            }
                        }
                    }
                }
                
                // next filter against the new candidates
                {
                    Iterator<HaplotypeCandidate> potentialCandidateIter =
                        potentialNewCandidateList.iterator();
                    while(potentialCandidateIter.hasNext())
                    {
                        HaplotypeCandidate currPotentialCandidate =
                            potentialCandidateIter.next();
                        
                        for(HaplotypeCandidate currNewCandidate: newHaplotypeCandidateList)
                        {
                            if(this.doesCandidate1ContainCandidate2(
                                    currNewCandidate,
                                    currPotentialCandidate))
                            {
                                // remove this candidate since it's contained
                                potentialCandidateIter.remove();
                                break;
                            }
                        }
                    }
                }
                
                // clear the haplotype candidate list... then add the new
                // stuff one-by-one
                haplotypeCandidateList.clear();
                
                // add the filtered "potential" candidates
                haplotypeCandidateList.addAll(potentialNewCandidateList);
                
                // copy the new list over to the candidate list
                haplotypeCandidateList.addAll(
                        newHaplotypeCandidateList);
            }
            
            // Now we just need to clean up all of the haplotype candidates
            for(HaplotypeCandidate currHaplotypeCandidate: haplotypeCandidateList)
            {
                int indexDifference =
                    snpCount - currHaplotypeCandidate.getStartingSnpIndex();
                if(indexDifference >= windowSizeInSnps)
                {
                    haplotypeList.add(IntervalScanningHaplotypeEstimatorTest.createHaplotypeBlockFromCandidate(
                            currHaplotypeCandidate,
                            snpCount,
                            sortedChromosomes));
                }
            }
        }
        else
        {
            // too few chromosomes to discover haplotypes...
            // just log it
            System.out.println(
                    "no haplotypes to find since we only have " +
                    chromosomes.size() + " chromosomes");
        }
        
        return haplotypeList;
    }
    
    private static PartitionedInterval createHaplotypeBlockFromCandidate(
            HaplotypeCandidate candidate,
            int terminationIndex,
            StrainChromosome[] chromosomes)
    {
        StrainChromosome firstChromosome =
            chromosomes[candidate.getStrainBitSet().nextSetBit(0)];
        int startingSnpIndex =
            candidate.getStartingSnpIndex();
        int endingSnpIndex = terminationIndex - 1;
        assert startingSnpIndex <= endingSnpIndex;
        
        SingleNucleotidePolymorphism startingSnp =
            firstChromosome.getSingleNucleotidePolymorphisms()[startingSnpIndex];
        SingleNucleotidePolymorphism endingSnp =
            firstChromosome.getSingleNucleotidePolymorphisms()[endingSnpIndex];
        long startingPositionBp =
            startingSnp.getPositionInBasePairs();
        long endingPositionBp =
            endingSnp.getPositionInBasePairs();
        assert startingPositionBp <= endingPositionBp;
        
        PartitionedInterval acceptedHaplotype = new PartitionedInterval(
                firstChromosome.getChromosomeNumber(),
                startingSnp.getPositionInBasePairs(),
                1 + endingPositionBp - startingPositionBp,
                candidate.getStrainBitSet());
        
        return acceptedHaplotype;
    }

    /**
     * Check to see whether or not the 1st haplotype contains the 2nd.
     * @param haplotype1
     *          the 1st haplotype
     * @param haplotype2
     *          the 2nd haplotype
     * @return
     *          true if the 1st candidate contains the second
     */
    private boolean doesCandidate1ContainCandidate2(
            HaplotypeCandidate haplotype1,
            HaplotypeCandidate haplotype2)
    {
        if(haplotype1.getStartingSnpIndex() <= haplotype2.getStartingSnpIndex())
        {
            return SetUtilities.isSubset(
                    haplotype2.getStrainBitSet(),
                    haplotype1.getStrainBitSet());
        }
        
        return false;
    }
    
    /**
     * Perform some basic validity checks on the chromosomes.
     * @param chromosomes
     *          the chromosomes to check
     * @throws IllegalArgumentException
     *          if the chromosomes are not correctly structured for this
     *          estimator
     */
    private static void validateChromosomes(StrainChromosome[] chromosomes)
    throws IllegalArgumentException
    {
        if(chromosomes.length >= 2)
        {
            int chromosomeNumber = chromosomes[0].getChromosomeNumber();
            int snpCount =
                chromosomes[0].getSingleNucleotidePolymorphisms().length;
            
            // make sure that all of the chromosome #'s  and SNP counts match up
            for(StrainChromosome currChromosome: chromosomes)
            {
                int currChromosomeNumber = currChromosome.getChromosomeNumber();
                if(currChromosomeNumber != chromosomeNumber)
                {
                    throw new IllegalArgumentException(
                            "Cannot estimate haplotypes for different " +
                            "chromosome numbers: " + chromosomeNumber  +
                            " and " + currChromosomeNumber);
                }
                
                int currSnpCount =
                    currChromosome.getSingleNucleotidePolymorphisms().length;
                if(currSnpCount != snpCount)
                {
                    throw new IllegalArgumentException(
                            "Cannot estimate haplotypes for chromosomeList " +
                            "with different SNP counts: " + snpCount +
                            " and " + currSnpCount);
                }
            }
            
            // make sure that all strains are unique
            for(int i = 0; i < chromosomes.length; i++)
            {
                for(int j = i + 1; j < chromosomes.length; j++)
                {
                    if(chromosomes[i].getStrainName().equals(chromosomes[j].getStrainName()))
                    {
                        throw new IllegalArgumentException(
                                "Cannot estimate haplotypes where the strain " +
                                "names match: " + chromosomes[i].getStrainName());
                    }
                }
            }
        }
    }
    
    /**
     * Group all of the chromosomes by their snp values at the given index.
     * If the chromosomes to group passed in are sorted, then the returned
     * groups should also be sorted
     * @param chromosomesToGroup
     *          the chromosomes that we're going to separate into different groups
     * @param snpIndexToGroupAt
     *          the index that we're grouping at
     * @return
     *          the groups
     */
    private static BitSet[] groupChromosomesBySnpTypeAt(
            StrainChromosome[] chromosomesToGroup,
            int snpIndexToGroupAt)
    {
        EnumMap<SnpType, BitSet> chromosomeGroupingsMap =
            new EnumMap<SnpType, BitSet>(SnpType.class);
        
        // put each chromosome in the right group
        for(int i = 0; i < chromosomesToGroup.length; i++)
        {
            StrainChromosome currChromosome = chromosomesToGroup[i];
            
            SnpType currSnpType =
                currChromosome.getSingleNucleotidePolymorphisms()[snpIndexToGroupAt].getSnpType();
            BitSet matchingGroup = chromosomeGroupingsMap.get(currSnpType);
            
            if(matchingGroup == null)
            {
                matchingGroup = new BitSet();
                chromosomeGroupingsMap.put(
                        currSnpType,
                        matchingGroup);
            }
            matchingGroup.set(i);
        }
        
        // convert to array
        return chromosomeGroupingsMap.values().toArray(
                new BitSet[chromosomeGroupingsMap.size()]);
    }
}
