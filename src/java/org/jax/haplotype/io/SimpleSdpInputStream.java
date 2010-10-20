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
 * A simple {@link SdpInputStream} that sits on a list of SNP input
 * streams
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SimpleSdpInputStream implements SdpInputStream
{
    private final SnpInputStream[] snpInputStreams;
    
    private final String[] strainNames;
    
    /**
     * Constructor
     * @param sdpStrainNames
     *          the strain names
     * @param snpInputStreams
     *          the SNP input streams that are the underlying data sources
     *          for this SDP input stream. The bit indices used from
     *          {@link #getNextSdp()} will match the stream ordering
     * @throws IOException
     *          if we get an exception during stream initialization
     */
    public SimpleSdpInputStream(
            String[] sdpStrainNames,
            SnpInputStream[] snpInputStreams)
    throws IOException
    {
        this.snpInputStreams = snpInputStreams;
        
        for(int i = 1; i < this.snpInputStreams.length; i++)
        {
            SnpInputStream prevSnpStream = this.snpInputStreams[i - 1];
            SnpInputStream currSnpStream = this.snpInputStreams[i];
            
            if(prevSnpStream.getSnpCount() != currSnpStream.getSnpCount())
            {
                throw new IllegalArgumentException(
                        "SNP counts from all of the streams should match!");
            }
        }
        
        this.strainNames = sdpStrainNames;
    }

    /**
     * {@inheritDoc}
     */
    public BitSet getNextSdp() throws IOException
    {
        BitSet nextSdp = new BitSet(this.snpInputStreams.length);
        
        for(int i = 0; i < this.snpInputStreams.length; i++)
        {
            if(this.snpInputStreams[i].getNextSnp())
            {
                nextSdp.set(i);
            }
        }
        
        return nextSdp;
    }

    /**
     * {@inheritDoc}
     */
    public long getSdpCount() throws IOException
    {
        return this.snpInputStreams[0].getSnpCount();
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSdp() throws IOException
    {
        return this.snpInputStreams[0].hasNextSnp();
    }
    
    /**
     * {@inheritDoc}
     */
    public String[] getSdpStrainNames() throws IOException
    {
        return this.strainNames;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.snpInputStreams[0].getReadDirection();
    }
}
