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
package org.jax.haplotype.phylogeny.inference;

import org.jax.geneticutil.data.SimpleBasePairInterval;

/**
 * A {@link SimpleBasePairInterval} that also has a {@link SdpInclusionHierarchy}
 * @author Keith Sheppard
 */
/*package-protected*/ class SdpInclusionInterval extends SimpleBasePairInterval
{
    /**
     * every {@link java.io.Serializable} should have one of these
     */
    private static final long serialVersionUID = 3109196877163974973L;
    
    private final SdpInclusionHierarchy sdpInclusionHierarchy;

    /**
     * Constructor
     * @param chromosomeNumber
     *          see {@link #getChromosomeNumber()}
     * @param startInBasePairs
     *          see {@link #getStartInBasePairs()}
     * @param extentInBasePairs
     *          see {@link #getExtentInBasePairs()}
     * @param sdpInclusionHierarchy
     *          see {@link #getSdpInclusionHierarchy()}
     */
    public SdpInclusionInterval(
            int chromosomeNumber,
            long startInBasePairs,
            long extentInBasePairs,
            SdpInclusionHierarchy sdpInclusionHierarchy)
    {
        super(chromosomeNumber, startInBasePairs, extentInBasePairs);
        this.sdpInclusionHierarchy = sdpInclusionHierarchy;
    }
    
    /**
     * Getter for the SDP inclusion hierarchy
     * @return the sdpInclusionHierarchy
     */
    public SdpInclusionHierarchy getSdpInclusionHierarchy()
    {
        return this.sdpInclusionHierarchy;
    }
}
