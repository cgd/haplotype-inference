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
