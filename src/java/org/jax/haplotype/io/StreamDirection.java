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
