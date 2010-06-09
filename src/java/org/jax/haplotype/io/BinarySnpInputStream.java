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
import java.io.InputStream;

/**
 * Reads bit packed snp data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class BinarySnpInputStream implements SnpInputStream
{
    private final InputStream inputStream;
    
    private final StreamDirection readDirection;
    
    private final long snpCount;
    
    private long index = 0;
    
    private int currByte = 0;

    /**
     * The byte masks used to pull out the bits we want
     */
    /*package protected*/ static final byte[] BYTE_MASKS = new byte[] {
        0x01,
        0x02,
        0x04,
        0x08,
        0x10,
        0x20,
        0x40,
        (byte)0x80};
    
    /**
     * Constructor
     * @param inputStream
     *          the input stream to use
     * @throws IOException
     *          if the snp count read fails
     */
    public BinarySnpInputStream(InputStream inputStream) throws IOException
    {
        this.inputStream = inputStream;
        this.readDirection = StreamDirection.byteToStreamDirection(
                (byte)inputStream.read());
        DataInputStream dataInputStream = new DataInputStream(inputStream);
        this.snpCount = dataInputStream.readLong();
    }
    
    /**
     * {@inheritDoc}
     */
    public boolean getNextSnp() throws IOException
    {
        if(this.index == this.snpCount)
        {
            throw new IOException();
        }
        
        int modIndex = (int)(this.index % 8);
        if(modIndex == 0)
        {
            this.currByte = this.inputStream.read();
        }
        
        int bitVal =
            this.currByte &
            BYTE_MASKS[modIndex];
        this.index++;
        
        return bitVal != 0;
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
    public boolean hasNextSnp() throws IOException
    {
        return this.index < this.snpCount;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.readDirection;
    }
}
