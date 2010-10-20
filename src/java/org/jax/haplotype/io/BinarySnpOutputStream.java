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

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/**
 * A bit-packing snp output stream
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class BinarySnpOutputStream implements SnpOutputStream
{
    /**
     * the stream to write to
     */
    private final OutputStream outputStream;
    
    private byte currentByte;
    
    private long index;

    private final long snpCount;

    private final StreamDirection writeDirection;
    
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
     * @param writeDirection
     *          the write direction that this stream uses
     * @param outputStream
     *          the stream to write the bits to
     * @param snpCount
     *          the snp count
     * @throws IOException
     *          if we fail to write initial bytes to stream
     */
    public BinarySnpOutputStream(
            StreamDirection writeDirection,
            OutputStream outputStream,
            long snpCount) throws IOException
    {
        this.writeDirection = writeDirection;
        this.outputStream = outputStream;
        this.snpCount = snpCount;
        
        this.outputStream.write(StreamDirection.streamDirectionToByte(
                writeDirection));
        DataOutputStream dataOutputStream = new DataOutputStream(outputStream);
        dataOutputStream.writeLong(snpCount);
    }

    /**
     * {@inheritDoc}
     */
    public void writeSnp(boolean snpMatchesReference) throws IOException
    {
        int modIndex = (int)(this.index % 8);
        if(snpMatchesReference)
        {
            this.currentByte |= BYTE_MASKS[modIndex];
        }
        
        if(modIndex == 7 || this.index == this.snpCount - 1)
        {
            this.commitByteToStream();
        }
        
        this.index++;
    }
    
    /**
     * Commit the byte to file
     * @throws IOException
     *          if the commit fails
     */
    private void commitByteToStream() throws IOException
    {
        this.outputStream.write(this.currentByte);
        this.currentByte = 0x00;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getWriteDirection()
    {
        return this.writeDirection;
    }
    
    /**
     * a little main test
     * @param args
     *          don't care
     */
    public static void main(String[] args)
    {
        System.out.println((byte)0xFF);
        System.out.println(0xFFFF & (byte)0xFF);
        System.out.println(0xFF);
        System.out.println(0xFFFF & 0xFF);
        System.out.println((0xFFFF & (byte)0xFF) & 0xFF);
    }
}
