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
package org.jax.haplotype.io;

import java.io.IOException;
import java.util.BitSet;

/**
 * An SDP stream that normalizes the SDP's by always ensuring that the
 * minority allele is set to 1, and if it's a 50/50 split then
 * the 1st strain is set to 0.
 * @author Keith Sheppard
 */
public class MinorityNormalizedSdpInputStream implements SdpInputStream
{
    private final SdpInputStream delegateSdpInputStream;
    
    private final int strainCount;
    
    /**
     * Constructor
     * @param delegateSdpInputStream
     *          the delegate stream that whose SDP output we'll normalize
     * @throws IOException
     *          if our delegate throws an {@link IOException} during
     *          initialization
     */
    public MinorityNormalizedSdpInputStream(
            SdpInputStream delegateSdpInputStream)
    throws IOException
    {
        this.delegateSdpInputStream = delegateSdpInputStream;
        this.strainCount =
            this.delegateSdpInputStream.getSdpStrainNames().length;
    }

    /**
     * {@inheritDoc}
     */
    public BitSet getNextSdp() throws IOException
    {
        BitSet nextSdp = this.delegateSdpInputStream.getNextSdp();
        
        // if there are more 1's than 0's OR if the 1's count is equal to
        // the 0's count AND the 1st bit is on, then we need to flip all
        // of the bits to normalize
        int doubleOnesCount = nextSdp.cardinality() * 2;
        if(doubleOnesCount > this.strainCount ||
           (doubleOnesCount == this.strainCount && nextSdp.get(0)))
        {
            nextSdp.flip(0, this.strainCount);
        }
        
        return nextSdp;
    }

    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.delegateSdpInputStream.getReadDirection();
    }

    /**
     * {@inheritDoc}
     */
    public long getSdpCount() throws IOException
    {
        return this.delegateSdpInputStream.getSdpCount();
    }

    /**
     * {@inheritDoc}
     */
    public String[] getSdpStrainNames() throws IOException
    {
        return this.delegateSdpInputStream.getSdpStrainNames();
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSdp() throws IOException
    {
        return this.delegateSdpInputStream.hasNextSdp();
    }
}
