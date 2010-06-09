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
import java.util.List;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.util.io.CharacterDelimitedParser;
import org.jax.util.io.IllegalFormatException;


/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommaSeparatedSnpBlockParser
{
    private static final String DEFAULT_CHROMOSOME_HEADER =
            "chromosome";
    
    private static final String DEFAULT_BLOCK_START_HEADER =
        "blockStartingPositionInBasePairs";
    
    private static final String DEFAULT_BLOCK_EXTENT_HEADER =
        "blockExtentInBasePairs";
    
    private final CharacterDelimitedParser commaSeparatedParser = new CharacterDelimitedParser();
    
    /**
     * Parse SNP blocks from the given input stream using default column
     * header names
     * @param is
     *          the input stream to parse
     * @return
     *          the SNP block list
     * @throws IOException
     *          if parsing the stream fails
     * @throws IllegalFormatException
     *          if the format is bad
     */
    public List<BasePairInterval> parseSnpBlocks(
            InputStream is) throws IOException, IllegalFormatException
    {
        return this.parseSnpBlocks(
                is,
                DEFAULT_CHROMOSOME_HEADER,
                DEFAULT_BLOCK_START_HEADER,
                DEFAULT_BLOCK_EXTENT_HEADER);
    }
    
    /**
     * Parse SNP blocks from the given input stream using the given column
     * header names
     * @param is
     *          the input stream to parse
     * @param chromosomeHeader
     *          the chromosome header name
     * @param blockStartHeader
     *          the block starting position column name
     * @param blockExtentHeader
     *          the block extent column name
     * @return
     *          the SNP block list
     * @throws IOException
     *          if parsing the stream fails
     * @throws IllegalFormatException
     *          if the format is bad
     */
    public List<BasePairInterval> parseSnpBlocks(
            InputStream is,
            String chromosomeHeader,
            String blockStartHeader,
            String blockExtentHeader) throws IOException, IllegalFormatException
    {
        BufferedReader bufferedReader = new BufferedReader(
                new InputStreamReader(is));
        
        String currLine = bufferedReader.readLine();
        while(currLine != null && this.commaSeparatedParser.isCommentOrEmpty(currLine))
        {
            currLine = bufferedReader.readLine();
        }
        if(currLine == null)
        {
            throw new IllegalFormatException("Can't parse empty file");
        }
        
        String[] currLineTokens;
        int chromosomeIndex;
        int blockStartIndex;
        int blockExtentIndex;
        if(chromosomeHeader == null)
        {
            if(blockStartHeader != null || blockExtentHeader != null)
            {
                throw new IllegalArgumentException(
                        "If chromosomeHeader is null the other headers " +
                        "should be null too");
            }
            
            chromosomeIndex = 0;
            blockStartIndex = 1;
            blockExtentIndex = 2;
        }
        else
        {
            currLineTokens = this.commaSeparatedParser.parseCharacterDelimitedLine(
                    currLine);
            
            chromosomeIndex = -1;
            blockStartIndex = -1;
            blockExtentIndex = -1;
            for(int i = 0; i < currLineTokens.length; i++)
            {
                if(currLineTokens[i].equalsIgnoreCase(chromosomeHeader))
                {
                    if(chromosomeIndex != -1)
                    {
                        throw new IllegalFormatException(
                                "Chromosome header found in both column " +
                                (chromosomeIndex + 1) + " and " +
                                (i + 1));
                    }
                    
                    chromosomeIndex = i;
                }
                else if(currLineTokens[i].equalsIgnoreCase(blockStartHeader))
                {
                    if(blockStartIndex != -1)
                    {
                        throw new IllegalFormatException(
                                "Block start header found in both column " +
                                (blockStartIndex + 1) + " and " +
                                (i + 1));
                    }
                    
                    blockStartIndex = i;
                }
                else if(currLineTokens[i].equalsIgnoreCase(blockExtentHeader))
                {
                    if(blockExtentIndex != -1)
                    {
                        throw new IllegalFormatException(
                                "Block extent header found in both column " +
                                (blockExtentIndex + 1) + " and " +
                                (i + 1));
                    }
                    
                    blockExtentIndex = i;
                }
            }
            
            if(chromosomeIndex == -1)
            {
                throw new IllegalFormatException(
                        "Failed to find " + chromosomeHeader);
            }
            if(blockStartIndex == -1)
            {
                throw new IllegalFormatException(
                        "Failed to find " + blockStartHeader);
            }
            if(blockExtentIndex == -1)
            {
                throw new IllegalFormatException(
                        "Failed to find " + blockExtentHeader);
            }
        }
        
        List<BasePairInterval> snpBlockList =
            new ArrayList<BasePairInterval>();
        while((currLineTokens = this.commaSeparatedParser.parseCharacterDelimitedLine(bufferedReader)) != null)
        {
            int chromosomeNumber = GenotypeParser.chromoNumStringToChromoNum(
                    currLineTokens[chromosomeIndex]);
            long blockStart = Long.parseLong(currLineTokens[blockStartIndex]);
            long blockExtent = Long.parseLong(currLineTokens[blockExtentIndex]);
            
            // convert a negative block into a real block
            if(blockExtent < 0)
            {
                blockExtent = -blockExtent;
                blockStart -= blockExtent;
            }
            
            snpBlockList.add(new SimpleBasePairInterval(
                    chromosomeNumber,
                    blockStart,
                    blockExtent));
        }
        
        return snpBlockList;
    }
}
