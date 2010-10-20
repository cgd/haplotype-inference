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
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IntersectionIntervalExpression implements SnpIntervalExpression
{
    private final SnpIntervalExpression expression1;
    
    private final SnpIntervalExpression expression2;
    
    /**
     * Create a new intersection expression
     * @param expression1
     *          the 1st expression to intersect
     * @param expression2
     *          the 2nd expression to intersect
     */
    public IntersectionIntervalExpression(
            SnpIntervalExpression expression1,
            SnpIntervalExpression expression2)
    {
        this.expression1 = expression1;
        this.expression2 = expression2;
    }
    
    /**
     * {@inheritDoc}
     */
    public BasePairInterval[] evaluateExpression()
    {
        BasePairInterval[] intervals1 = this.expression1.evaluateExpression();
        BasePairInterval[] intervals2 = this.expression2.evaluateExpression();
        
        return IntersectionIntervalExpression.intersect(intervals1, intervals2);
    }
    
    /**
     * Intersect intersect the given snp intervals
     * @param intervals1
     *          the 1st set of intervals
     * @param intervals2
     *          the 2nd set of intervals
     * @return
     *          the result of the intersection
     */
    public static BasePairInterval[] intersect(
            BasePairInterval[] intervals1,
            BasePairInterval[] intervals2)
    {
        ArrayList<BasePairInterval> intersectionIntervals =
            new ArrayList<BasePairInterval>(
                    Math.min(intervals1.length, intervals2.length));
        
        int intervals1Cursor = 0;
        int intervals2Cursor = 0;
        
        while(intervals1Cursor < intervals1.length && intervals2Cursor < intervals2.length)
        {
            BasePairInterval currInterval1 = intervals1[intervals1Cursor];
            BasePairInterval currInterval2 = intervals2[intervals2Cursor];
            BasePairInterval currIntersection =
                IntersectionIntervalExpression.intersect(currInterval1, currInterval2);
            
            if(currIntersection != null)
            {
                intersectionIntervals.add(currIntersection);
            }
            
            if(currInterval1.getEndInBasePairs() < currInterval2.getEndInBasePairs())
            {
                intervals1Cursor++;
            }
            else
            {
                intervals2Cursor++;
            }
        }
        
        return intersectionIntervals.toArray(
                new BasePairInterval[intersectionIntervals.size()]);
    }

    /**
     * Perform intersection on the two given snp intervals
     * @param interval1
     *          the 1st snp interval
     * @param interval2
     *          the 2nd snp interval
     * @return
     *          the intersection
     */
    private static BasePairInterval intersect(
            BasePairInterval interval1,
            BasePairInterval interval2)
    {
        long maxStart = Math.max(
                interval1.getStartInBasePairs(),
                interval2.getStartInBasePairs());
        long minEnd = Math.min(
                interval1.getEndInBasePairs(),
                interval2.getEndInBasePairs());
        
        int chrom1 = interval1.getChromosomeNumber();
        int chrom2 = interval2.getChromosomeNumber();
        if(maxStart <= minEnd && chrom1 == chrom2)
        {
            return new SimpleBasePairInterval(
                    chrom1,
                    maxStart,
                    1 + (minEnd - maxStart));
        }
        else
        {
            return null;
        }
    }
}
