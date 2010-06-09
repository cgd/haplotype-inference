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

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class UnionIntervalExpression implements SnpIntervalExpression
{
    private final SnpIntervalExpression expression1;
    
    private final SnpIntervalExpression expression2;
    
    /**
     * Create a new union expression that does a union operation on the
     * given snp intervals
     * @param expression1
     *          the 1st expression to perform a union operation on
     * @param expression2
     *          the 2nd expression to perform a union operation on
     */
    public UnionIntervalExpression(
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
        
        return UnionIntervalExpression.calculateUnion(
                intervals1,
                intervals2);
    }
    
    /**
     * calculate the union of the two given interval arrays
     * @param intervals1
     *          the 1st interval list
     * @param intervals2
     *          the 2nd interval list
     * @return
     *          the union of the intervals
     */
    public static BasePairInterval[] calculateUnion(
            BasePairInterval[] intervals1,
            BasePairInterval[] intervals2)
    {
        ArrayList<BasePairInterval> unionIntervals =
            new ArrayList<BasePairInterval>(
                    Math.min(intervals1.length, intervals2.length));
        
        int intervals1Cursor = 0;
        int intervals2Cursor = 0;
        
        while(intervals1Cursor < intervals1.length && intervals2Cursor < intervals2.length)
        {
            int chromosome1 = intervals1[intervals1Cursor].getChromosomeNumber();
            int chromosome2 = intervals2[intervals2Cursor].getChromosomeNumber();
            
            if(chromosome1 == chromosome2)
            {
                long unionStart = Math.min(
                        intervals1[intervals1Cursor].getStartInBasePairs(),
                        intervals2[intervals2Cursor].getStartInBasePairs());
                long unionEnd = -1;
                boolean foundGap = false;
                
                // iterate through until we hit a gap
                while(!foundGap && intervals1Cursor < intervals1.length && intervals2Cursor < intervals2.length)
                {
                    if(UnionIntervalExpression.intervalsAreTouchingOrOverlapping(
                            intervals1[intervals1Cursor],
                            intervals2[intervals2Cursor]))
                    {
                        if(intervals1[intervals1Cursor].getEndInBasePairs() < intervals2[intervals2Cursor].getEndInBasePairs())
                        {
                            intervals1Cursor++;
                        }
                        else
                        {
                            intervals2Cursor++;
                        }
                    }
                    else
                    {
                        // we found a gap
                        foundGap = true;
                        if(intervals1[intervals1Cursor].getEndInBasePairs() < intervals2[intervals2Cursor].getEndInBasePairs())
                        {
                            unionEnd = intervals1[intervals1Cursor].getEndInBasePairs();
                            intervals1Cursor++;
                        }
                        else
                        {
                            unionEnd = intervals2[intervals2Cursor].getEndInBasePairs();
                            intervals2Cursor++;
                        }
                    }
                }
                
                if(!foundGap)
                {
                    // we ran into the end of one of the intervals before we hit a gap
                    if(intervals1Cursor == intervals1.length)
                    {
                        unionEnd = Math.max(
                                intervals1[intervals1.length - 1].getEndInBasePairs(),
                                intervals2[intervals2Cursor].getEndInBasePairs());
                        intervals2Cursor++;
                    }
                    else
                    {
                        unionEnd = Math.max(
                                intervals1[intervals1Cursor].getEndInBasePairs(),
                                intervals2[intervals2.length - 1].getEndInBasePairs());
                        intervals1Cursor++;
                    }
                }
                
                unionIntervals.add(new SimpleBasePairInterval(
                        chromosome1,
                        unionStart,
                        1 + unionEnd - unionStart));
            }
        }
        
        // clean up
        for(; intervals1Cursor < intervals1.length; intervals1Cursor++)
        {
            unionIntervals.add(intervals1[intervals1Cursor]);
        }
        for(; intervals2Cursor < intervals2.length; intervals2Cursor++)
        {
            unionIntervals.add(intervals2[intervals2Cursor]);
        }
        
        return unionIntervals.toArray(
                new SimpleBasePairInterval[unionIntervals.size()]);
    }
    
    private static boolean intervalsAreTouchingOrOverlapping(
            BasePairInterval interval1,
            BasePairInterval interval2)
    {
        long maxStart = Math.max(
                interval1.getStartInBasePairs(),
                interval2.getStartInBasePairs());
        long minEnd = Math.min(
                interval1.getEndInBasePairs(),
                interval2.getEndInBasePairs());
        
        return maxStart - 1 <= minEnd;
    }
}
