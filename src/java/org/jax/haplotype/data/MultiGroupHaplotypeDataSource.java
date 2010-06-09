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
