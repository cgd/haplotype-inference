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

import java.util.Random;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.SnpType;
import org.jax.geneticutil.data.StrainChromosome;
import org.junit.Assert;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ExpressionTestUtil
{
    public static BasePairInterval[] createRandomBlocks(int blockCount)
    {
        Random rand = new Random();
        BasePairInterval[] randomBlocks = new BasePairInterval[blockCount];
        
        long lastTail = 0;
        for(int i = 0; i < randomBlocks.length; i++)
        {
            randomBlocks[i] = new SimpleBasePairInterval(
                    lastTail + rand.nextInt(10),
                    rand.nextInt(30) + 1);
            lastTail = randomBlocks[i].getEndInBasePairs() + 2;
        }
        
        return randomBlocks;
    }
    
    public static void validateBlocks(BasePairInterval[] blocks)
    {
        for(BasePairInterval currBlock: blocks)
        {
            Assert.assertTrue(currBlock.getStartInBasePairs() >= 0);
            Assert.assertTrue(currBlock.getExtentInBasePairs() >= 1);
            Assert.assertTrue(
                currBlock.getEndInBasePairs() ==
                currBlock.getStartInBasePairs() + currBlock.getExtentInBasePairs() - 1);
        }
        
        for(int i = 1; i < blocks.length; i++)
        {
            BasePairInterval currBlock = blocks[i];
            BasePairInterval prevBlock = blocks[i - 1];
            
            Assert.assertTrue(
                    prevBlock.getEndInBasePairs() < currBlock.getStartInBasePairs() - 1);
        }
    }
    
    public static StrainChromosome[] generateRandomChromosomes(
            int numChromosomes,
            int snpsPerChromosome,
            int chromosomeNumber,
            Random random)
    {
        StrainChromosome[] chromosomes = new StrainChromosome[numChromosomes];
        
        for(int chromosomeIndex = 0; chromosomeIndex < chromosomes.length; chromosomeIndex++)
        {
            chromosomes[chromosomeIndex] = new StrainChromosome(
                    Integer.toString(chromosomeIndex),
                    chromosomeNumber);
            chromosomes[chromosomeIndex].setSingleNucleotidePolymorphisms(
                    new SingleNucleotidePolymorphism[snpsPerChromosome]);
        }
        
        for(int snpIndex = 0; snpIndex < snpsPerChromosome; snpIndex++)
        {
            int rand1 =
                random.nextInt(SnpType.values().length);
            SnpType snpType1 =
                SnpType.values()[rand1];
            
            int rand2 =
                random.nextInt(SnpType.values().length - 1);
            if(rand2 >= rand1)
            {
                rand2++;
            }
            SnpType snpType2 =
                SnpType.values()[rand2];
            
            for(StrainChromosome currChromosome: chromosomes)
            {
                SingleNucleotidePolymorphism currSnp =
                    new SingleNucleotidePolymorphism(
                            random.nextBoolean() ? snpType1 : snpType2,
                            snpIndex);
                SingleNucleotidePolymorphism[] currSnps =
                    currChromosome.getSingleNucleotidePolymorphisms();
                currSnps[snpIndex] = currSnp;
            }
        }
        
        return chromosomes;
    }

    /**
     * @param blocks
     * @param snpIndex
     * @return
     */
    public static boolean blocksContainSnp(
            BasePairInterval[] blocks, int snpIndex)
    {
        for(BasePairInterval currBlock: blocks)
        {
            if(currBlock.getStartInBasePairs() <= snpIndex && currBlock.getEndInBasePairs() >= snpIndex)
            {
                return true;
            }
        }
        
        return false;
    }
}
