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

import java.io.DataInputStream;
import java.io.IOException;

/**
 * Stream for reading SNP positions stored in a binary format
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class BinarySnpPositionInputStream implements SnpPositionInputStream
{
    private final DataInputStream inputStream;
    
    private final StreamDirection readDirection;
    
    private final long snpCount;
    
    private final long extentInBasePairs;
    
    private final long minSnpPositionInBasePairs;
    
    private final int chromosomeNumber;
    
    private long index = 0L;

    /**
     * Constructor
     * @param inputStream
     *          the input stream to use
     * @throws IOException
     *          if the snp count read fails
     */
    public BinarySnpPositionInputStream(DataInputStream inputStream)
    throws IOException
    {
        this.inputStream = inputStream;
        this.readDirection = StreamDirection.byteToStreamDirection(
                (byte)inputStream.read());
        this.chromosomeNumber = this.inputStream.readInt();
        this.minSnpPositionInBasePairs = this.inputStream.readLong();
        this.extentInBasePairs = this.inputStream.readLong();
        this.snpCount = this.inputStream.readLong();
    }
    
    /**
     * {@inheritDoc}
     */
    public long getNextSnpPositionInBasePairs() throws IOException
    {
        this.index++;
        return this.inputStream.readLong();
    }
    
    /**
     * {@inheritDoc}
     */
    public long getSnpCount() throws IOException
    {
        return this.snpCount;
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean hasNextSnpPosition() throws IOException
    {
        return this.index < this.snpCount;
    }
    
    /**
     * {@inheritDoc}
     */
    public long getExtentInBasePairs()
    {
        return this.extentInBasePairs;
    }
    
    /**
     * {@inheritDoc}
     */
    public long getStartInBasePairs()
    {
        return this.minSnpPositionInBasePairs;
    }
    
    /**
     * {@inheritDoc}
     */
    public int getChromosomeNumber()
    {
        return this.chromosomeNumber;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.readDirection;
    }
}
