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
 * An SDP input stream whose 1/0 values are normalized vs a reference
 * SNP stream
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ReferenceNormalizedSdpInputStream implements SdpInputStream
{
    private final SnpInputStream referenceSnpStream;
    
    private final SdpInputStream comparisonSdpStream;
    
    private final int strainCount;
    
    /**
     * Constructor
     * @param referenceSnpStream
     *          the reference snp stream
     * @param comparisonSdpStream
     *          the comparison SDP stream
     * @throws IOException
     *          if we get an IO exception during initialization
     */
    public ReferenceNormalizedSdpInputStream(
            SnpInputStream referenceSnpStream,
            SdpInputStream comparisonSdpStream)
    throws IOException
    {
        this.referenceSnpStream = referenceSnpStream;
        this.comparisonSdpStream = comparisonSdpStream;
        
        this.strainCount = this.comparisonSdpStream.getSdpStrainNames().length;
    }

    /**
     * {@inheritDoc}
     */
    public BitSet getNextSdp() throws IOException
    {
        BitSet nextSdp = this.comparisonSdpStream.getNextSdp();
        if(!this.referenceSnpStream.getNextSnp())
        {
            nextSdp.flip(0, this.strainCount);
        }
        
        return nextSdp;
    }

    /**
     * {@inheritDoc}
     */
    public long getSdpCount() throws IOException
    {
        return this.comparisonSdpStream.getSdpCount();
    }

    /**
     * {@inheritDoc}
     */
    public String[] getSdpStrainNames() throws IOException
    {
        return this.comparisonSdpStream.getSdpStrainNames();
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSdp() throws IOException
    {
        return this.comparisonSdpStream.hasNextSdp();
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.referenceSnpStream.getReadDirection();
    }
}
