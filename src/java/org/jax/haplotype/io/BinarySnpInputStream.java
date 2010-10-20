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
