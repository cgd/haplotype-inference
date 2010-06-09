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
