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
