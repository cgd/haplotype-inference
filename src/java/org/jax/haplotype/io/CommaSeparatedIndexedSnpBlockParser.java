/*
 * Copyright (c) 2008 The Jackson Laboratory
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.jax.haplotype.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.jax.geneticutil.data.IndexedSnpInterval;
import org.jax.util.io.CharacterDelimitedParser;

/**
 * This parser reads in indexed SNP blocks. The expected format is:<br/>
 * <code>
 * chromosomeName,blockStartIndex,blockEndIndex
 * </code>
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommaSeparatedIndexedSnpBlockParser
{
    private final CharacterDelimitedParser commaSeparatedParser =
        new CharacterDelimitedParser(',', '#');
    
    /**
     * Read in snp blocks from the given input stream
     * @param inputStream
     *          the input stream to read snp blocks from
     * @return
     *          the SNP blocks
     * @throws IOException
     *          if there is some problem reading data from the stream
     */
    public List<IndexedSnpInterval> parseSnpBlocksFromStream(
            InputStream inputStream)
            throws IOException
    {
        BufferedReader br = new BufferedReader(new InputStreamReader(
                inputStream));
        List<IndexedSnpInterval> snpBlocks = new ArrayList<IndexedSnpInterval>();
        String[] currLine;
        while((currLine = this.commaSeparatedParser.parseCharacterDelimitedLine(br)) != null)
        {
            if(currLine.length != 3)
            {
                throw new IOException(
                        "Expected 3 comma separated values, but read " +
                        currLine.length + " in " + Arrays.toString(currLine));
            }
            else
            {
                IndexedSnpInterval currBlock = new IndexedSnpInterval(
                        Integer.parseInt(currLine[1]),
                        Integer.parseInt(currLine[2]));
                snpBlocks.add(currBlock);
            }
        }
        
        return snpBlocks;
    }
}
