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
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.List;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.util.TextWrapper;

/**
 * A class for writing SNP blocks to a comma-separated data stream
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommaSeparatedSnpBlockWriter
{
    private static final String COMMENT_PREFIX = "# ";
    
    private static final int HEADER_COMMENT_WRAP_MARGIN =
        120 - COMMENT_PREFIX.length();
    
    /**
     * Write the given SNP blocks to the given output stream and tack
     * on the given header comment
     * @param outputStream
     *          the output stream
     * @param snpBlocks
     *          the SNP blocks
     * @param headerComment
     *          the header comment. this comment does not have to have a
     *          '#' prefix since this function will add one
     * @throws IOException
     *          if we catch an {@link java.io.IOException} from the
     *          given output stream
     */
    public void writeSnpBlocksToStream(
            OutputStream outputStream,
            List<BasePairInterval> snpBlocks,
            String headerComment)
    throws IOException
    {
        PrintStream printStream = new PrintStream(
                outputStream);
        
        // write out the header comment
        if(headerComment != null)
        {
            String[] wrappedHeaderComment = TextWrapper.wrapText(
                    headerComment,
                    HEADER_COMMENT_WRAP_MARGIN);
            for(String headerCommentLine: wrappedHeaderComment)
            {
                headerCommentLine = COMMENT_PREFIX + headerCommentLine;
                printStream.println(headerCommentLine);
            }
            
            // write out the table header
            printStream.println(
                    "chromosome,intervalStartingPositionInBasePairs," +
                    "intervalExtentInBasePairs");
        }
        
        for(BasePairInterval currBlock: snpBlocks)
        {
            printStream.println(
                    currBlock.getChromosomeNumber() + "," +
                    currBlock.getStartInBasePairs() + "," +
                    currBlock.getExtentInBasePairs());
        }
        printStream.flush();
    }
}
