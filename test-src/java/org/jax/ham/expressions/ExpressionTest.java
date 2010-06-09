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

package org.jax.ham.expressions;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Random;
import java.util.Set;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.expressions.FunctionalExpressionParser;
import org.jax.haplotype.expressions.IntersectionIntervalExpression;
import org.jax.haplotype.expressions.MalformedExpressionException;
import org.jax.haplotype.expressions.NotSnpPositionEvaluator;
import org.jax.haplotype.expressions.SnpPositionExpression;
import org.jax.haplotype.expressions.SnpPositionsMatchEvaluator;
import org.jax.haplotype.expressions.UnionIntervalExpression;
import org.jax.haplotype.io.GenotypeParser;
import org.junit.Assert;
import org.junit.Test;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ExpressionTest
{
    @Test
    public void testIntersectionAndUnion()
    {
        BasePairInterval[] intersection;
        BasePairInterval[] block1;
        BasePairInterval[] block2;
        
        block1 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        0,
                        100)};
        block2 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        100,
                        100)};
        intersection = IntersectionIntervalExpression.intersect(
                block1,
                block2);
        Assert.assertTrue(intersection.length == 0);
        
        block1 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        0,
                        101)};
        block2 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        100,
                        100)};
        intersection = IntersectionIntervalExpression.intersect(
                block1,
                block2);
        Assert.assertTrue(intersection.length == 1);
        Assert.assertTrue(intersection[0].getExtentInBasePairs() == 1);
        Assert.assertTrue(intersection[0].getStartInBasePairs() == 100);
        
        block1 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        0,
                        100),
                new SimpleBasePairInterval(
                        101,
                        3),
                new SimpleBasePairInterval(
                        105,
                        1000)};
        block2 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        100,
                        200)};
        intersection = IntersectionIntervalExpression.intersect(
                block1,
                block2);
        System.out.println("intersection length: " + intersection.length);
        Assert.assertTrue(intersection.length == 2);
        Assert.assertEquals(
                intersection[0],
                new SimpleBasePairInterval(
                        101,
                        3));
        System.out.println("intersection[1]: " + intersection[1].toString());
        Assert.assertEquals(
                intersection[1],
                new SimpleBasePairInterval(
                        105,
                        195));
        
        block1 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        0,
                        100),
                new SimpleBasePairInterval(
                        101,
                        3),
                new SimpleBasePairInterval(
                        105,
                        1000)};
        block2 = new BasePairInterval[] {
                new SimpleBasePairInterval(
                        2,
                        1103)};
        intersection = IntersectionIntervalExpression.intersect(
                block1,
                block2);
        Assert.assertTrue(intersection.length == 3);
        
        Random random = new Random();
        for(int i = 0; i < 10000; i++)
        {
            int randNumSnps;
            if(random.nextBoolean())
            {
                randNumSnps = random.nextInt(20) + 1;
            }
            else
            {
                randNumSnps = random.nextInt(2000) + 1;
            }
            StrainChromosome[] randomChromosomes1;
            {
                int randNumChromo = random.nextInt(3) + 2;
                int randChromoNum = random.nextInt(50);
                
                randomChromosomes1 = ExpressionTestUtil.generateRandomChromosomes(
                        randNumChromo,
                        randNumSnps,
                        randChromoNum,
                        random);
            }
            
            StrainChromosome[] randomChromosomes2;
            {
                int randNumChromo = random.nextInt(3) + 2;
                int randChromoNum = random.nextInt(50);
                
                randomChromosomes2 = ExpressionTestUtil.generateRandomChromosomes(
                        randNumChromo,
                        randNumSnps,
                        randChromoNum,
                        random);
            }
            
            SnpPositionsMatchEvaluator positionsMatchEvaluator =
                new SnpPositionsMatchEvaluator();
            SnpPositionExpression matchExpression =
                new SnpPositionExpression(
                        randomChromosomes1,
                        positionsMatchEvaluator);
            NotSnpPositionEvaluator positionsDontMatchEvaluator =
                new NotSnpPositionEvaluator(new SnpPositionsMatchEvaluator());
            SnpPositionExpression notMatchExpression =
                new SnpPositionExpression(
                        randomChromosomes2,
                        positionsDontMatchEvaluator);
            
            BasePairInterval[] notMatchBlocks =
                notMatchExpression.evaluateExpression();
            BasePairInterval[] matchBlocks =
                matchExpression.evaluateExpression();
            ExpressionTestUtil.validateBlocks(notMatchBlocks);
            ExpressionTestUtil.validateBlocks(matchBlocks);
            
            System.out.println("==========");
            System.out.println("notMatchLen=" + notMatchBlocks.length);
            System.out.println("matchLen=" + matchBlocks.length);
            
            IntersectionIntervalExpression intersectionExpression =
                new IntersectionIntervalExpression(
                        matchExpression,
                        notMatchExpression);
            BasePairInterval[] intersectionBlocks =
                intersectionExpression.evaluateExpression();
            ExpressionTestUtil.validateBlocks(intersectionBlocks);
            System.out.println("intersectionLen=" + intersectionBlocks.length);
            
            
            UnionIntervalExpression unionExpression =
                new UnionIntervalExpression(
                        notMatchExpression,
                        matchExpression);
            BasePairInterval[] unionBlocks =
                unionExpression.evaluateExpression();
            ExpressionTestUtil.validateBlocks(unionBlocks);
            System.out.println("unionLen=" + unionBlocks.length);
            
            for(int j = 0; j < notMatchBlocks.length; j++)
            {
                BasePairInterval currNotMatchBlock =
                    notMatchBlocks[j];
                for(int k = 0; k < currNotMatchBlock.getExtentInBasePairs(); k++)
                {
                    // TODO reimplement test
//                    boolean matchesContain =
//                        ExpressionTestUtil.blocksContainSnp(
//                                matchBlocks,
//                                currNotMatchBlock.getStartInBasePairs() + k);
//                    boolean intersectionContains =
//                        ExpressionTestUtil.blocksContainSnp(
//                                intersectionBlocks,
//                                currNotMatchBlock.getStartInBasePairs() + k);
//                    boolean unionContains =
//                        ExpressionTestUtil.blocksContainSnp(
//                                unionBlocks,
//                                currNotMatchBlock.getStartInBasePairs() + k);
//                    Assert.assertTrue(unionContains);
//                    if(matchesContain)
//                    {
//                        Assert.assertTrue(intersectionContains);
//                    }
//                    else
//                    {
//                        Assert.assertFalse(intersectionContains);
//                    }
                }
            }
            
            for(int j = 0; j < matchBlocks.length; j++)
            {
                BasePairInterval currMatchBlock =
                    matchBlocks[j];
                for(int k = 0; k < currMatchBlock.getExtentInBasePairs(); k++)
                {
                    // TODO reimplement test
//                    boolean notMatchesContain =
//                        ExpressionTestUtil.blocksContainSnp(
//                                notMatchBlocks,
//                                currMatchBlock.getStartInBasePairs() + k);
//                    boolean intersectionContains =
//                        ExpressionTestUtil.blocksContainSnp(
//                                intersectionBlocks,
//                                currMatchBlock.getStartInBasePairs() + k);
//                    boolean unionContains =
//                        ExpressionTestUtil.blocksContainSnp(
//                                unionBlocks,
//                                currMatchBlock.getStartInBasePairs() + k);
//                    Assert.assertTrue(unionContains);
//                    if(notMatchesContain)
//                    {
//                        Assert.assertTrue(intersectionContains);
//                    }
//                    else
//                    {
//                        Assert.assertFalse(intersectionContains);
//                    }
                }
            }
        }
    }
    
    @Test
    public void testDataAndExpressionParsing() throws FileNotFoundException, IOException, MalformedExpressionException
    {
        GenotypeParser parser = new GenotypeParser();
        Set<StrainChromosome> chromosomes =
//            parser.parseGenotypeFromStream(new FileInputStream("test-data.csv"));
            parser.parseGenotypeFromStream(
                    ExpressionTest.class.getResourceAsStream("/chromosome_random_1.csv"));
        SnpPositionsMatchEvaluator matchEvaluator = new SnpPositionsMatchEvaluator();
        NotSnpPositionEvaluator notMatchEvaluator = new NotSnpPositionEvaluator(matchEvaluator);
        SnpPositionExpression matchOne = new SnpPositionExpression(
                this.getChromosomesNamed(chromosomes, new String[] {"C57BL/6J", "C3H/HeJ"}),
                matchEvaluator);
        SnpPositionExpression matchTwo = new SnpPositionExpression(
                this.getChromosomesNamed(chromosomes, new String[] {"SEG/Pas", "SJL/J", "SM/J"}),
                notMatchEvaluator);
        SnpPositionExpression matchThree = new SnpPositionExpression(
                this.getChromosomesNamed(chromosomes, new String[] {"NZB/BlNJ", "NZW/LacJ"}),
                matchEvaluator);
        
        UnionIntervalExpression union = new UnionIntervalExpression(matchTwo, matchThree);
        IntersectionIntervalExpression intersection = new IntersectionIntervalExpression(matchOne, union);
        BasePairInterval[] manualResults =
            intersection.evaluateExpression();
        System.out.println("manual snps: " + manualResults.length);
        
        FunctionalExpressionParser expParser = new FunctionalExpressionParser(chromosomes);
        String expressionString =
            "intersection(C57BL/6J=C3H/HeJ, union(not(SEG/Pas=SJL/J=SM/J), NZB/BlNJ=NZW/LacJ))";
        BasePairInterval[] parsedResults =
            expParser.parseExpression(expressionString).evaluateExpression();
        System.out.println("parsed snps: " + parsedResults.length);
        
        Assert.assertTrue(Arrays.equals(manualResults, parsedResults));
    }

    /**
     * @param strings
     * @return
     */
    private StrainChromosome[] getChromosomesNamed(Set<StrainChromosome> chromosomeSet, String[] names)
    {
        StrainChromosome[] chromosomes = new StrainChromosome[names.length];
        OUTER_LOOP:
        for(int i = 0; i < names.length; i++)
        {
            for(StrainChromosome chromosome: chromosomeSet)
            {
                if(chromosome.getStrainName().equals(names[i]))
                {
                    chromosomes[i] = chromosome;
                    continue OUTER_LOOP;
                }
            }
            System.out.println("Coult not find: " + names[i]);
            System.out.println("chromosomes:");
            for(StrainChromosome chromosome: chromosomeSet)
            {
                System.out.print(" " + chromosome.getStrainName());
            }
            Assert.fail();
        }
        
        return chromosomes;
    }
}
