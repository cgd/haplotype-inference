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
