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

/**
 * An input stream for reading SNP positions
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface SnpPositionInputStream
{
    /**
     * Get the chromosome number for this position stream
     * @return
     *          the chromosome number
     */
    public int getChromosomeNumber();
    
    /**
     * Get the lowest base pair value
     * @return
     *          the lowest base pair value
     */
    public long getStartInBasePairs();
    
    /**
     * Get the extent base-pair value
     * @return
     *          the extent base-pair value
     */
    public long getExtentInBasePairs();
    
    /**
     * Get the SNP count
     * @return
     *          the snp count
     * @throws IOException
     *          if we fail to get the snp count
     */
    public long getSnpCount() throws IOException;
    
    /**
     * Move to the next snp and get its position in base pairs
     * @return
     *          the next snp position
     * @throws IOException
     *          if we get an exception reading data
     */
    public long getNextSnpPositionInBasePairs() throws IOException;
    
    /**
     * Determine if there are any more snps left
     * @return
     *          true if there are any more snps left
     * @throws IOException
     *          if the read fails for some reason
     */
    public boolean hasNextSnpPosition() throws IOException;
    
    /**
     * Get the chromosome read direction used by this stream
     * @return
     *          the read direction
     * @throws IOException
     *          if we fail to get the read direction
     */
    public StreamDirection getReadDirection() throws IOException;
}
