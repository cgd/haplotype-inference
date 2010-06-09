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

package org.jax.haplotype.inference;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.geneticutil.data.PartitionedIntervalSet;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.geneticutil.data.SimplePartitionedIntervalSet;

/**
 * A class for extracting equivalence classes from a list of haplotype blocks
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HaplotypeEquivalenceClassCreator
{
    /**
     * Convert the list of haplotype blocks into a list of equivalence classes
     * by finding all of the unique SDP (Strain Distribution Pattern)
     * @param haplotypeBlocks
     *          the haplotype blocks to convert
     * @return
     *          the equivalence classes
     */
    public static List<PartitionedIntervalSet> createEquivalenceClassesFromBlocks(
            List<PartitionedInterval> haplotypeBlocks)
    {
        Map<BitSet, List<PartitionedInterval>> haplotypeEquivMap =
            new HashMap<BitSet, List<PartitionedInterval>>();
        for(PartitionedInterval currHaplotypeBlock: haplotypeBlocks)
        {
            List<PartitionedInterval> currHaplotypeList =
                haplotypeEquivMap.get(currHaplotypeBlock.getStrainBitSet());
            if(currHaplotypeList == null)
            {
                currHaplotypeList = new ArrayList<PartitionedInterval>();
                haplotypeEquivMap.put(
                        currHaplotypeBlock.getStrainBitSet(),
                        currHaplotypeList);
            }
            
            currHaplotypeList.add(currHaplotypeBlock);
        }
        
        List<PartitionedIntervalSet> equivalenceClassList =
            new ArrayList<PartitionedIntervalSet>(haplotypeEquivMap.size());
        Iterator<Entry<BitSet, List<PartitionedInterval>>> equivalenceEntryIter =
            haplotypeEquivMap.entrySet().iterator();
        while(equivalenceEntryIter.hasNext())
        {
            Entry<BitSet, List<PartitionedInterval>> currEntry =
                equivalenceEntryIter.next();
            List<PartitionedInterval> currHapList = currEntry.getValue();
            
            BasePairInterval[] newSnpBlocks =
                new BasePairInterval[currHapList.size()];
            for(int blockIndex = 0; blockIndex < newSnpBlocks.length; blockIndex++)
            {
                PartitionedInterval currHaplotypeBlock = currHapList.get(blockIndex);
                newSnpBlocks[blockIndex] = new SimpleBasePairInterval(
                        currHaplotypeBlock.getChromosomeNumber(),
                        currHaplotypeBlock.getStartInBasePairs(),
                        currHaplotypeBlock.getExtentInBasePairs());
            }
            
            equivalenceClassList.add(new SimplePartitionedIntervalSet(
                    currEntry.getKey(),
                    newSnpBlocks));
            
            equivalenceEntryIter.remove();
        }
        
        return equivalenceClassList;
    }
}
