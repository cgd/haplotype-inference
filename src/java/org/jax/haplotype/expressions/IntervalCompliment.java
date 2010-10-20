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

package org.jax.haplotype.expressions;

import java.util.ArrayList;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;

/**
 * Class for calculating the compliment of a set of intervals
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntervalCompliment
{
    /**
     * Calculate the complement of the input interval
     * @param intervalsToCompliment
     *          the interval that we're getting a compliment for. these
     *          should be in order
     * @param chromosomeNumber
     *          the chromosome number to compliment. all given intervals should
     *          be for this chromosome
     * @param startingPositionBp
     *          the starting position in base pairs
     * @param totalExtentBp
     *          the total extent in base pairs
     * @return
     *          the compliment of the input... given the full range values
     */
    public static BasePairInterval[] calculateCompliment(
            BasePairInterval[] intervalsToCompliment,
            int chromosomeNumber,
            long startingPositionBp,
            long totalExtentBp)
    {
        if(intervalsToCompliment.length == 0)
        {
            BasePairInterval interval = new SimpleBasePairInterval(
                    chromosomeNumber,
                    startingPositionBp,
                    totalExtentBp);
            return new BasePairInterval[] {interval};
        }
        else
        {
            long endingPositionBp = startingPositionBp + totalExtentBp - 1;
            
            ArrayList<BasePairInterval> complimentIntervals = new ArrayList<BasePairInterval>(
                    intervalsToCompliment.length + 1);
            
            // the 1st interval needs special treatment
            BasePairInterval firstInterval = intervalsToCompliment[0];
            if(startingPositionBp > firstInterval.getStartInBasePairs())
            {
                long endBp = firstInterval.getStartInBasePairs() - 1;
                complimentIntervals.add(new SimpleBasePairInterval(
                        chromosomeNumber,
                        startingPositionBp,
                        1 + endBp - startingPositionBp));
            }
            else if(startingPositionBp > firstInterval.getStartInBasePairs())
            {
                throw new IllegalArgumentException(
                        "intervals to compliment starts before the given " +
                        "starting position");
            }
            
            // deal with the middle intervals
            for(int i = 0; i < intervalsToCompliment.length - 1; i++)
            {
                BasePairInterval interval1 = intervalsToCompliment[i];
                BasePairInterval interval2 = intervalsToCompliment[i + 1];
                
                long startBp = interval1.getEndInBasePairs() + 1;
                long endBp = interval2.getStartInBasePairs() - 1;
                if(startBp > endBp)
                {
                    throw new IllegalArgumentException(
                            "overlap detected at index: " + i);
                }
                
                complimentIntervals.add(new SimpleBasePairInterval(
                        chromosomeNumber,
                        startBp,
                        1 + endBp - startBp));
            }
            
            // the last interval needs special treatment
            BasePairInterval lastInterval =
                intervalsToCompliment[intervalsToCompliment.length - 1];
            if(endingPositionBp > lastInterval.getEndInBasePairs())
            {
                long startBp = lastInterval.getEndInBasePairs() + 1;
                complimentIntervals.add(new SimpleBasePairInterval(
                        chromosomeNumber,
                        startBp,
                        1 + endingPositionBp - startBp));
            }
            else if(endingPositionBp < lastInterval.getEndInBasePairs())
            {
                throw new IllegalArgumentException(
                        "intervals to compliment ends after the given " +
                        "ending position");
            }
            
            // ok... all done
            return complimentIntervals.toArray(
                    new BasePairInterval[complimentIntervals.size()]);
        }
    }
}
