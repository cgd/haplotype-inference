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
