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

package org.jax.haplotype.expressions;

import java.util.ArrayList;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.geneticutil.data.StrainChromosome;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SnpPositionExpression implements SnpIntervalExpression
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            SnpPositionExpression.class.getName());
    
    private final StrainChromosome[] chromosomesToCompare;

    private final SnpPositionEvaluator snpPositionEvaluator;
    
    /**
     * Constructor
     * @param chromosomesToCompare
     *          the chromosomes we're comparing
     * @param snpPositionEvaluator
     *          the evaluator to use
     */
    public SnpPositionExpression(
            StrainChromosome[] chromosomesToCompare,
            SnpPositionEvaluator snpPositionEvaluator)
    {
        this.chromosomesToCompare = chromosomesToCompare;
        this.snpPositionEvaluator = snpPositionEvaluator;
    }

    /**
     * {@inheritDoc}
     */
    public BasePairInterval[] evaluateExpression()
    {
        if(this.chromosomesToCompare == null)
        {
            throw new NullPointerException("chromosomes are null!");
        }
        else if(this.chromosomesToCompare.length < 2)
        {
            throw new IllegalStateException(
                    "no chromosomes to compare! length = " +
                    this.chromosomesToCompare.length);
        }
        else
        {
            ArrayList<BasePairInterval> snpIntervals =
                new ArrayList<BasePairInterval>();
            
            int minSnpCount = this.chromosomesToCompare[0].getSingleNucleotidePolymorphisms().length;
            {
                int chromosomeNumber = this.chromosomesToCompare[0].getChromosomeNumber();
                for(int i = 1; i < this.chromosomesToCompare.length; i++)
                {
                    StrainChromosome currChromosome = this.chromosomesToCompare[i];
                    if(chromosomeNumber != currChromosome.getChromosomeNumber())
                    {
                        LOG.warning("chromosome numbers don't agree... continuing");
                    }
                    
                    if(minSnpCount != currChromosome.getSingleNucleotidePolymorphisms().length)
                    {
                        LOG.warning("SNP counts don't agree on chromosomes... continuing");
                        minSnpCount = Math.min(
                                minSnpCount,
                                currChromosome.getSingleNucleotidePolymorphisms().length);
                    }
                }
                
                int currIntervalStartIndex = -1;
                for(int currSnpIndex = 0; currSnpIndex < minSnpCount; currSnpIndex++)
                {
                    if(this.snpPositionEvaluator.evaluate(this.chromosomesToCompare, currSnpIndex))
                    {
                        if(currIntervalStartIndex == -1)
                        {
                            currIntervalStartIndex = currSnpIndex;
                        }
                    }
                    else
                    {
                        if(currIntervalStartIndex != -1)
                        {
                            snpIntervals.add(new SimpleBasePairInterval(
                                    currIntervalStartIndex,
                                    currSnpIndex - currIntervalStartIndex));
                            currIntervalStartIndex = -1;
                        }
                    }
                }
                
                if(currIntervalStartIndex != -1)
                {
                    snpIntervals.add(new SimpleBasePairInterval(
                            currIntervalStartIndex,
                            minSnpCount - currIntervalStartIndex));
                }
            }
            
            // convert the intervals to an array and return it
            return snpIntervals.toArray(
                    new SimpleBasePairInterval[snpIntervals.size()]);
        }
    }
}
