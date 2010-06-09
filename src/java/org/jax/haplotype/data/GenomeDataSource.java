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

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Map;
import java.util.Set;

import org.jax.util.datastructure.SequenceUtilities;

/**
 * A collection of chromosome data sources
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenomeDataSource implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 5960172600154047538L;
    
    private final Map<Integer, ? extends ChromosomeDataSource> chromosomeDataSources;
    private final String name;
    private final String ncbiBuildVersion;
    
    /**
     * Constructor
     * @param name
     *          the name of this genome data source
     * @param ncbiBuildVersion
     *          the NCBI build version
     * @param chromosomeDataSources
     *          the chromosome data sources
     */
    public GenomeDataSource(
            String name,
            String ncbiBuildVersion,
            Map<Integer, ? extends ChromosomeDataSource> chromosomeDataSources)
    {
        this.name = name;
        this.ncbiBuildVersion = ncbiBuildVersion;
        this.chromosomeDataSources = chromosomeDataSources;
    }
    
    /**
     * Getter for the NCBI build version
     * @return
     *          the build version string
     */
    public String getNcbiBuildVersion()
    {
        return this.ncbiBuildVersion;
    }
    
    /**
     * Getter for the available strains
     * @return
     *          the available strains
     */
    public Set<String> getAvailableStrains()
    {
        if(this.chromosomeDataSources.isEmpty())
        {
            return Collections.emptySet();
        }
        else
        {
            // we're assuming here that all of the chromosomes have the same
            // set of strains
            ChromosomeDataSource anyChromosomeDataSource =
                this.chromosomeDataSources.values().iterator().next();
            return anyChromosomeDataSource.getAvailableStrains();
        }
    }
    
    /**
     * Getter for the available chromosomes
     * @return
     *          the available chromosomes
     */
    public int[] getAvailableChromosomes()
    {
        int[] availableChromosomes = SequenceUtilities.toIntArray(
                new ArrayList<Integer>(this.chromosomeDataSources.keySet()));
        Arrays.sort(availableChromosomes);
        return availableChromosomes;
    }
    
    /**
     * Get the chromosomes mapped by chromosome number
     * @return the chromosome data sources
     */
    public Map<Integer, ? extends ChromosomeDataSource> getChromosomeDataSources()
    {
        return this.chromosomeDataSources;
    }
    
    /**
     * Getter for the name
     * @return the name
     */
    public String getName()
    {
        return this.name;
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
