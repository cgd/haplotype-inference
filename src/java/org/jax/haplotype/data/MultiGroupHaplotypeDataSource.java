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

package org.jax.haplotype.data;

import java.io.Serializable;
import java.util.List;
import java.util.Set;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.geneticutil.data.PartitionedInterval;

/**
 * This is a haplotype data source that returns
 * {@link MultiPartitionedInterval}s rather than {@link PartitionedInterval}s
 * (effectively you have multiple groups per partition so that you can do
 * something like an F-test rather than just a t-test).
 * 
 * TODO I think the name I chose for this class is confusing and should be
 *      changed to something that is more intuitive
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface MultiGroupHaplotypeDataSource extends Serializable
{
    /**
     * Getter for the name of this haplotype data source
     * @return
     *          the name of the haplotype data source (can be null)
     */
    public String getName();
    
    /**
     * Get the haplotype data
     * @param chromosome
     *          the chromosome
     * @param strainsToAccept
     *          the strains to allow through
     * @return
     *          the list of haplotype blocks given the strain and chromosome
     *          filters
     */
    public List<MultiPartitionedInterval> getHaplotypeData(
            int chromosome,
            Set<String> strainsToAccept);
    
    /**
     * Get the available strains in the data
     * @return
     *          the strains
     */
    public Set<String> getAvailableStrains();
    
    /**
     * Getter for the available chromosome numbers
     * @return
     *          the chromosomes
     */
    public int[] getAvailableChromosomes();
}
