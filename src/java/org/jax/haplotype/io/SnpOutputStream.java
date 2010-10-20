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
 * Simple interface for writing snp data
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface SnpOutputStream
{
    /**
     * Write a boolean value to stream saying that the current snp matches
     * the reference snp
     * @param snpMatchesReference
     *          value that determines whether or not we match the reference
     * @throws IOException
     *          if we have trouble writing the data
     */
    public void writeSnp(boolean snpMatchesReference)
    throws IOException;
    
    /**
     * Get the direction that this stream writes in
     * @return
     *          this stream's write direction
     */
    public StreamDirection getWriteDirection();
}
