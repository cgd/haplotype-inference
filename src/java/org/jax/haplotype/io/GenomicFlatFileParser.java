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

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jax.util.io.CharacterDelimitedParser;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenomicFlatFileParser implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -3982158591427126252L;
    
    // Example Header 1: snp.ID,ChrID,build.36.bp.Position,Source,129S1/SvImJ,129S4/SvJae,129X1/SvJ
    // Example Header 2: snp ID,ChrID,build 36 bp Position,Source,129S1/SvImJ,129S4/SvJae,129X1/SvJ
    private static final Pattern CHROMOSOME_ID_HEADER_PATTERN = Pattern.compile(
            "^chr.*",
            Pattern.CASE_INSENSITIVE);
    private static final String CONFLICT_HEADER = "idx.validated.conflict.snp";
    private static final Pattern BASE_PAIR_HEADER_PATTERN = Pattern.compile(
        ".*position.*",
        Pattern.CASE_INSENSITIVE);
    
    /**
     * for matching #'s
     */
    private static final Pattern GREEDY_DIGIT_PATTERN = Pattern.compile("[0-9]+");
    
    /**
     * the number of annotation columns for this format (this is the number
     * of columns before we hit strain snp data)
     */
    private static final int DEFAULT_NUM_ANNOTATION_COLUMNS = 4;
    
    /**
     * our inner parser
     */
    protected final CharacterDelimitedParser characterDelimitedParser;
    
    /**
     * the number of annotation columns
     */
    protected final int numAnnotationColumns;
    
    /**
     * the column index of the chromosome ID
     */
    protected final int indexOfChromosomeId;
    
    /**
     * the column index for the base pair position
     */
    protected final int indexOfBasePairPosition;
    
    /**
     * Constructor
     */
    public GenomicFlatFileParser()
    {
        this(new CharacterDelimitedParser(), DEFAULT_NUM_ANNOTATION_COLUMNS);
    }
    
    /**
     * Constructor
     * @param characterDelimitedParser
     *          the parser
     * @param numAnnotationColumns
     *          the number of annotation columns
     */
    public GenomicFlatFileParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns)
    {
        this(characterDelimitedParser,
             numAnnotationColumns,
             -1,
             -1);
    }
    
    /**
     * Constructor
     * @param numAnnotationColumns
     *          the number of annotation columns
     * @param indexOfChromosomeId
     *          the chromosome ID index (-1 means don't know and we'll take a
     *          guess)
     * @param indexOfBasePairPosition
     *          the base pair position index (-1 means don't know and we'll
     *          take a guess)
     */
    public GenomicFlatFileParser(
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition)
    {
        this(new CharacterDelimitedParser(),
             numAnnotationColumns,
             indexOfChromosomeId,
             indexOfBasePairPosition);
    }
    
    /**
     * Constructor
     * @param characterDelimitedParser
     *          the parser
     * @param numAnnotationColumns
     *          the number of annotation columns
     * @param indexOfChromosomeId
     *          the chromosome ID index (-1 means don't know and we'll take a
     *          guess)
     * @param indexOfBasePairPosition
     *          the base pair position index (-1 means don't know and we'll
     *          take a guess)
     */
    public GenomicFlatFileParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition)
    {
        this.characterDelimitedParser = characterDelimitedParser;
        this.numAnnotationColumns = numAnnotationColumns;
        this.indexOfChromosomeId = indexOfChromosomeId;
        this.indexOfBasePairPosition = indexOfBasePairPosition;
    }
    
    /**
     * Parse the reader into a set of chromosomes.
     * @param bufferedReader
     *          the reader that we're parsing
     * @return
     *          the header info
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public HeaderInfo parseHeaderInfoFromReader(
            BufferedReader bufferedReader) throws IOException
    {
        // 1st we need the header
        String[] headerStrings = this.characterDelimitedParser.parseCharacterDelimitedLine(bufferedReader);
        if(headerStrings.length <= this.numAnnotationColumns)
        {
            throw new IOException(
                    "comma separated file contains fewer than " +
                    (this.numAnnotationColumns + 1) +
                    " columns, meaning that we have no data");
        }
        
        final int indexOfChromosomeId;
        if(this.indexOfChromosomeId == -1)
        {
            indexOfChromosomeId = guessIndexOfChromosomeIdFromHeader(
                    headerStrings,
                    this.numAnnotationColumns);
            if(indexOfChromosomeId < 0)
            {
                throw new IOException(
                        "could not find chromosome ID within the first " +
                        this.numAnnotationColumns + " columns");
            }
        }
        else
        {
            indexOfChromosomeId = this.indexOfChromosomeId;
        }
        
        final int indexOfBasePairPosition;
        if(this.indexOfBasePairPosition == -1)
        {
            indexOfBasePairPosition = guessIndexOfBasePairPositionFromHeader(
                    headerStrings,
                    this.numAnnotationColumns);
            if(indexOfBasePairPosition < 0)
            {
                throw new IOException(
                        "could not find base pair header within the first " +
                        this.numAnnotationColumns + " columns");
            }
        }
        else
        {
            indexOfBasePairPosition = this.indexOfBasePairPosition;
        }
        
        // i'm assuming that a conflict column has to be the last column
        boolean hasConflictColumn = CONFLICT_HEADER.equalsIgnoreCase(
                headerStrings[headerStrings.length - 1]);
        int rowLengthWithoutConflictColumn =
            hasConflictColumn ? headerStrings.length - 1 : headerStrings.length;
        
        return new HeaderInfo(
                headerStrings,
                indexOfChromosomeId,
                indexOfBasePairPosition,
                rowLengthWithoutConflictColumn);
    }
    
    /**
     * Convert the given chromosome number string to a number.
     * @param chromoNumString
     *          the string we're converting to a number The only thing
     *          we require is that there's a number somewhere in the string, or
     *          that it ends in X, Y or M. If it ends in any of the letters,
     *          the mouse chromosome numbering is assumed (20, 21, 22)
     * @return
     *          the number
     */
    public static int chromoNumStringToChromoNum(String chromoNumString)
    {
        chromoNumString = chromoNumString.trim();
        
        Matcher chromoNumMatcher = GREEDY_DIGIT_PATTERN.matcher(chromoNumString);
        if(chromoNumMatcher.find())
        {
            String numOnlyString = chromoNumMatcher.group();
            return Integer.parseInt(numOnlyString);
        }
        else
        {
            if(chromoNumString.endsWith("X") || chromoNumString.endsWith("x"))
            {
                return 20;
            }
            else if(chromoNumString.endsWith("Y") || chromoNumString.endsWith("y"))
            {
                return 21;
            }
            else if(chromoNumString.endsWith("M") || chromoNumString.endsWith("m"))
            {
                return 22;
            }
            else
            {
                throw new IllegalArgumentException(
                        "could not parse chromosome number string \"" +
                        chromoNumString + "\"");
            }
        }
    }
    
    /**
     * Get the index of the base pair position from the header
     * @param header
     *          the header
     * @return
     *          the index
     */
    public static int guessIndexOfBasePairPositionFromHeader(
            String[] header)
    {
        return guessIndexOfBasePairPositionFromHeader(header, header.length);
    }
    
    /**
     * Get the index of the base pair position from the header
     * @param header
     *          the header
     * @param numberOfColumnsToCheck
     *          the number of columns that we should search through before
     *          giving up
     * @return
     *          the index
     */
    private static int guessIndexOfBasePairPositionFromHeader(
            String[] header,
            int numberOfColumnsToCheck)
    {
        // try to find the index of the BP position
        for(int i = 0; i < header.length && i < numberOfColumnsToCheck; i++)
        {
            Matcher currMatcher = BASE_PAIR_HEADER_PATTERN.matcher(
                    header[i]);
            if(currMatcher.matches())
            {
                // found it!
                return i;
            }
        }
        
        // guess we couldn't find it
        return -1;
    }

    /**
     * Get the index of the chromosome id column from the given header strings
     * @param header
     *          the array of header strings
     * @return
     *          the index of the chromosome ID column or -1 if it isn't found
     */
    public static int guessIndexOfChromosomeIdFromHeader(String[] header)
    {
        return guessIndexOfChromosomeIdFromHeader(header, header.length);
    }
    
    /**
     * Get the index of the chromosome id column from the given header strings
     * @param header
     *          the array of header strings
     * @param numberOfColumnsToCheck
     *          the number of columns that we should search through before
     *          giving up
     * @return
     *          the index of the chromosome ID column or -1 if it isn't found
     */
    private static int guessIndexOfChromosomeIdFromHeader(
            String[] header,
            int numberOfColumnsToCheck)
    {
        // try to find the index of the chromosome
        for(int i = 0; i < header.length && i < numberOfColumnsToCheck; i++)
        {
            Matcher currMatcher = CHROMOSOME_ID_HEADER_PATTERN.matcher(
                    header[i]);
            if(currMatcher.matches())
            {
                // found it!
                return i;
            }
        }
        
        // guess we couldn't find it
        return -1;
    }

    /**
     * Determines which strains exist in the data presented by the given
     * input stream
     * @param inputStream
     *          the input stream
     * @return
     *          the set of strain names
     * @throws IOException
     *          if we run into trouble trying to parse data from the given
     *          stream
     */
    public Set<String> parseAvailableStrainsFromStream(
            InputStream inputStream) throws IOException
    {
        return new HashSet<String>(this.parseAvailableStrainsInOrderFromStream(
                inputStream));
    }
    
    /**
     * Determines which strains exist in the data presented by the given
     * input stream
     * @param inputStream
     *          the input stream
     * @return
     *          the ordered list of strain names
     * @throws IOException
     *          if we run into trouble trying to parse data from the given
     *          stream
     */
    public List<String> parseAvailableStrainsInOrderFromStream(
            InputStream inputStream) throws IOException
    {
        BufferedReader bufferedReader = new BufferedReader(
                new InputStreamReader(inputStream));
        
        // 1st we need the header
        String[] headerStrings = this.characterDelimitedParser.parseCharacterDelimitedLine(bufferedReader);
        if(headerStrings.length <= this.numAnnotationColumns)
        {
            throw new IOException(
                    "comma separated file contains fewer than " +
                    (this.numAnnotationColumns + 1) +
                    " columns, meaning that we have no data");
        }
        
        // i'm assuming that a conflict column has to be the last column
        boolean hasConflictColumn = CONFLICT_HEADER.equalsIgnoreCase(headerStrings[headerStrings.length - 1]);
        
        // OK now build the strain set
        List<String> strainList = new ArrayList<String>();
        String[] currRow = new String[headerStrings.length];
        int rowLengthWithoutConflictColumn =
            hasConflictColumn ? currRow.length - 1 : currRow.length;
        for(int i = this.numAnnotationColumns; i < rowLengthWithoutConflictColumn; i++)
        {
            strainList.add(headerStrings[i]);
        }
        
        return strainList;
    }
    
    /**
     * Contains info about the header
     */
    public static final class HeaderInfo
    {
        private final String[] headerStrings;
        
        private final int chromosomeIdIndex;
        
        private final int basePairPositionIndex;
        
        private final int rowLengthWithoutConflictColumn;
        
        /**
         * Constructor
         * @param headerStrings
         *          the header strings
         * @param chromosomeIdIndex
         *          the chromosome id index
         * @param basePairPositionIndex
         *          the base pair position index
         * @param rowLengthWithoutConflictColumn
         *          the row length not counting a conflict column
         */
        public HeaderInfo(
                String[] headerStrings,
                int chromosomeIdIndex,
                int basePairPositionIndex,
                int rowLengthWithoutConflictColumn)
        {
            this.headerStrings = headerStrings;
            this.chromosomeIdIndex = chromosomeIdIndex;
            this.basePairPositionIndex = basePairPositionIndex;
            this.rowLengthWithoutConflictColumn = rowLengthWithoutConflictColumn;
        }

        /**
         * Getter for the header strings
         * @return
         *          the header strings
         */
        public String[] getHeaderStrings()
        {
            return this.headerStrings;
        }

        /**
         * Getter for the chromosome id column index
         * @return
         *          the chromosome id index
         */
        public int getChromosomeIdIndex()
        {
            return this.chromosomeIdIndex;
        }

        /**
         * Getter for the base pair position column index
         * @return
         *          the index
         */
        public int getBasePairPositionIndex()
        {
            return this.basePairPositionIndex;
        }

        /**
         * Getter for the row length not counting any conflict column (if
         * it exists)
         * @return
         *          the length
         */
        public int getRowLengthWithoutConflictColumn()
        {
            return this.rowLengthWithoutConflictColumn;
        }
    }
    
    /**
     * Extract the chromosome number from the give row
     * @param headerInfo
     *          the header information to use
     * @param row
     *          the row
     * @return
     *          the chromosome number
     */
    protected int getChromosomeNumber(
            HeaderInfo headerInfo,
            String[] row)
    {
        int chromosomeNum = GenotypeParser.chromoNumStringToChromoNum(
                row[headerInfo.getChromosomeIdIndex()]);
        return chromosomeNum;
    }
    
    /**
     * Extract the base pair position from the given row
     * @param headerInfo
     *          the header info to use
     * @param row
     *          the row
     * @return
     *          the position
     */
    protected long getBasePairPosition(
            HeaderInfo headerInfo,
            String[] row)
    {
        long currBasePairPosition = Long.parseLong(
                row[headerInfo.getBasePairPositionIndex()]);
        return currBasePairPosition;
    }
}
