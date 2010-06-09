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

package org.jax.haplotype.data;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.haplotype.inference.SlidingWindowMultiHaplotypeEstimator;

/**
 * A test class that uses the sliding window algorithm
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SlidingWindowHaplotypeDataSource implements MultiGroupHaplotypeDataSource
{
    /**
     * every {@link java.io.Serializable}
     */
    private static final long serialVersionUID = -2641977364858670897L;
    
    private final String name;
    private final GenomeDataSource genomeDataSource;
    private final SlidingWindowMultiHaplotypeEstimator slidingWindowEstimator;
    
    /**
     * Constructor
     * @param name
     *          the name to use for this test (should be unique among tests)
     * @param genomeDataSource
     *          the genome data source to use
     * @param windowSizeInSnps
     *          the window size to use in SNPs
     * @param stepWindowBySnps
     *          should we step the window forward one snp at a time? see
     *          {@link SlidingWindowMultiHaplotypeEstimator#getStepWindowBySnps()}
     */
    public SlidingWindowHaplotypeDataSource(
            String name,
            GenomeDataSource genomeDataSource,
            int windowSizeInSnps,
            boolean stepWindowBySnps)
    {
        this.name = name;
        this.slidingWindowEstimator = new SlidingWindowMultiHaplotypeEstimator(
                windowSizeInSnps,
                stepWindowBySnps);
        this.genomeDataSource = genomeDataSource;
    }
    
    /**
     * Getter for the name of this test
     * @return the name
     */
    public String getName()
    {
        return this.name;
    }
    
    /**
     * Getter for the genome data source
     * @return the genome data source
     */
    public GenomeDataSource getGenomeDataSource()
    {
        return this.genomeDataSource;
    }
    
    /**
     * {@inheritDoc}
     */
    public int[] getAvailableChromosomes()
    {
        return this.getGenomeDataSource().getAvailableChromosomes();
    }
    
    /**
     * {@inheritDoc}
     */
    public Set<String> getAvailableStrains()
    {
        return this.genomeDataSource.getAvailableStrains();
    }
    
    /**
     * {@inheritDoc}
     */
    public List<MultiPartitionedInterval> getHaplotypeData(
            int chromosome,
            Set<String> strainsToAcceptFilter)
    {
        // find the common set of strains and sort them
        if(strainsToAcceptFilter == null)
        {
            strainsToAcceptFilter = this.getAvailableStrains();
        }
        else
        {
            strainsToAcceptFilter = new HashSet<String>(strainsToAcceptFilter);
            strainsToAcceptFilter.retainAll(this.getAvailableStrains());
        }
        
        String[] strainsToAcceptArray = strainsToAcceptFilter.toArray(
                new String[strainsToAcceptFilter.size()]);
        Arrays.sort(strainsToAcceptArray);
        
        ChromosomeDataSource chromoDataSource =
            this.genomeDataSource.getChromosomeDataSources().get(
                    chromosome);
        
        try
        {
            List<MultiPartitionedInterval> haplotypeDataList =
                this.slidingWindowEstimator.estimateMultiHaplotypeBlocks(
                        chromoDataSource.getSdpInputStream(strainsToAcceptArray),
                        chromoDataSource.getSnpPositionInputStream());
            return haplotypeDataList;
        }
        catch(RuntimeException ex)
        {
            // allow runtime exceptions to pass through
            throw ex;
        }
        catch(Exception ex)
        {
            // this may be a bad thing to do. we're repackaging the exception
            // and throwing it as a runtime
            throw new RuntimeException(ex);
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return this.getName();
    }
}
