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

import org.jax.geneticutil.data.CompositeRealValuedBasePairInterval;


/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyTestResult
extends CompositeRealValuedBasePairInterval
implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 6855647783889139607L;

    private final PhylogenyInterval phylogenyInterval;
    
    /**
     * Constructor
     * @param phylogenyInterval
     *          the phylogeny interval that the p values apply to
     * @param pValue
     *          the p-values
     */
    public PhylogenyTestResult(
            PhylogenyInterval phylogenyInterval,
            double pValue)
    {
        super(phylogenyInterval.getInterval(), pValue);
        this.phylogenyInterval = phylogenyInterval;
    }
    
    /**
     * Getter for the phylogeny interval
     * @return the phylogenyInterval
     */
    public PhylogenyInterval getPhylogenyInterval()
    {
        return this.phylogenyInterval;
    }
    
    /**
     * Getter for the p-value
     * @return the pValue
     */
    public double getPValue()
    {
        return this.getRealValue();
    }
}
