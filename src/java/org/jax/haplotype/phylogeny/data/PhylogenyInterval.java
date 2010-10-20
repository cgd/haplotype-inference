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

import org.jax.geneticutil.data.BasePairInterval;

/**
 * An aggregation of {@link PhylogenyTreeNode} and
 * {@link BasePairInterval}
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyInterval implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 3699248954123833123L;

    private final PhylogenyTreeNode phylogeny;
    
    private final BasePairInterval interval;

    /**
     * Constructor
     * @param phylogeny
     *          the phylogeny tree node
     * @param interval
     *          the SNP interval that this phylogeny applies to
     */
    public PhylogenyInterval(
            PhylogenyTreeNode phylogeny,
            BasePairInterval interval)
    {
        this.phylogeny = phylogeny;
        this.interval = interval;
    }
    
    /**
     * Getter for the phylogeny tree
     * @return the phylogeny
     */
    public PhylogenyTreeNode getPhylogeny()
    {
        return this.phylogeny;
    }
    
    /**
     * Getter for the interval that this tree applies to
     * @return the interval
     */
    public BasePairInterval getInterval()
    {
        return this.interval;
    }
}
