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

/**
 * Indicates the direction we're streaming the genome/chromosome in
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public enum StreamDirection
{
    /**
     * stream from start to  end
     */
    FORWARD,
    
    /**
     * stream from end to start
     */
    REVERSE;
    
    /**
     * Convert the given byte value to a stream direction
     * @param streamDirectionByte
     *          the byte value
     * @return
     *          the byte value
     */
    public static StreamDirection byteToStreamDirection(final byte streamDirectionByte)
    {
        switch(streamDirectionByte)
        {
            case 0:
            {
                return FORWARD;
            }
            
            case 1:
            {
                return REVERSE;
            }
            
            default:
            {
                throw new IllegalArgumentException(
                        "Invalid StreamDirection byte value: " +
                        streamDirectionByte);
            }
        }
    }
    
    /**
     * Convert the given {@link StreamDirection} enum to a byte value
     * @param streamDirection
     *          the stream direction to convert
     * @return
     *          the matching byte value
     */
    public static byte streamDirectionToByte(final StreamDirection streamDirection)
    {
        switch(streamDirection)
        {
            case FORWARD:
            {
                return 0;
            }
            
            case REVERSE:
            {
                return 1;
            }
            
            default:
            {
                throw new IllegalArgumentException(
                        "don't know how to convert given stream direction " +
                        "to a byte value: " + streamDirection);
            }
        }
    }
}
