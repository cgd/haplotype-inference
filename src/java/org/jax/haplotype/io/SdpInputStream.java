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
 * An input stream for reading SNP data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface SdpInputStream
{
    /**
     * Getter for the SDP count
     * @return
     *          the sdp count
     * @throws IOException
     *          if we fail to read the next sdp
     */
    public long getSdpCount() throws IOException;
    
    /**
     * Getter for the next SDP
     * @return
     *          the next SDP
     * @throws IOException
     *          if the read fails
     */
    public BitSet getNextSdp() throws IOException;
    
    /**
     * Determine if there are any more SDPs
     * @return
     *          true if there are more SDPs
     * @throws IOException
     *          if the read fails
     */
    public boolean hasNextSdp() throws IOException;
    
    /**
     * Get the strain names for the SDPs. the ordering used in this array
     * should match up with the bit ordering used for {@link #getNextSdp()}
     * @return
     *          the strain names
     * @throws IOException
     *          if theres a read error
     */
    public String[] getSdpStrainNames() throws IOException;
    
    /**
     * Get the chromosome read direction used by this stream
     * @return
     *          the read direction
     * @throws IOException
     *          if we fail to get the read direction
     */
    public StreamDirection getReadDirection() throws IOException;
}
