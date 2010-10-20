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

package org.jax.haplotype.phylogeny.data;

import java.io.Serializable;

import org.jax.geneticutil.data.IndexedSnpInterval;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IndexedPhylogenyInterval implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 6181753762556468516L;

    private PhylogenyTreeNode phylogeny;
    
    private IndexedSnpInterval indexedInterval;
    
    /**
     * Constructor
     * @param indexedInterval
     *          see {@link #getIndexedInterval()}
     * @param phylogeny
     *          see {@link #getPhylogeny()}
     */
    public IndexedPhylogenyInterval(
            PhylogenyTreeNode phylogeny,
            IndexedSnpInterval indexedInterval)
    {
        this.indexedInterval = indexedInterval;
        this.phylogeny = phylogeny;
    }
    
    /**
     * Getter for the phylogeny
     * @return the phylogeny
     */
    public PhylogenyTreeNode getPhylogeny()
    {
        return this.phylogeny;
    }
    
    /**
     * Getter for the indexed interval
     * @return the indexed interval
     */
    public IndexedSnpInterval getIndexedInterval()
    {
        return this.indexedInterval;
    }
}
