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

import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.util.io.CharacterDelimitedParser;

/**
 * A parser that extracts data from the HMM state output
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HiddenMarkovModelStateParser extends GenomicFlatFileParser
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -6116214569515301017L;
    
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            HiddenMarkovModelStateParser.class.getName());
    
    private final boolean flattenMatchingHaplotypes;
    
    /**
     * Constructor
     */
    public HiddenMarkovModelStateParser()
    {
        super();
        
        this.flattenMatchingHaplotypes = true;
    }
    
    /**
     * Constructor
     * @param characterDelimitedParser
     *          the parser
     * @param numAnnotationColumns
     *          the number of annotation columns
     */
    public HiddenMarkovModelStateParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns)
    {
        super(characterDelimitedParser, numAnnotationColumns);
        
        this.flattenMatchingHaplotypes = true;
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
    public HiddenMarkovModelStateParser(
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition)
    {
        super(new CharacterDelimitedParser(),
             numAnnotationColumns,
             indexOfChromosomeId,
             indexOfBasePairPosition);
        
        this.flattenMatchingHaplotypes = true;
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
     * @param flattenMatchingHaplotypes
     *          should contiguous haplotypes be flattened into a single interval
     *          when their haplotype groups match exactly (this is usually a
     *          good thing since you lose no information and are left with a
     *          more compact representation)
     */
    public HiddenMarkovModelStateParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition,
            boolean flattenMatchingHaplotypes)
    {
        super(characterDelimitedParser,
              numAnnotationColumns,
              indexOfChromosomeId,
              indexOfBasePairPosition);
        
        this.flattenMatchingHaplotypes = flattenMatchingHaplotypes;
    }
    
    /**
     * Parse the reader into a map of chromosomes with haplotypes.
     * <BR/><BR/>
     * <pre>
     * <em>
     * NOTE: the strain indices used in the {@link MultiPartitionedInterval}s
     *       are determined by the sort order of the strain names in the files
     *       header row NOT the column ordering in the file. 
     * </em>
     * </pre>
     * @param bufferedReader
     *          the reader that we're parsing
     * @param strainsToParse
     *          the names of the strains we should parse... if null parse
     *          everything
     * @return
     *          the map of chromosomes with haplotypes
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Map<Integer, List<MultiPartitionedInterval>> parseHMMStatesFromReader(
            BufferedReader bufferedReader,
            Set<String> strainsToParse) throws IOException
    {
        HeaderInfo headerInfo = this.parseHeaderInfoFromReader(bufferedReader);
        
        // figure out what the sorted strain ordering is
        ArrayList<String> unsortedStrains = new ArrayList<String>();
        for(int i = this.numAnnotationColumns;
            i < headerInfo.getRowLengthWithoutConflictColumn();
            i++)
        {
            unsortedStrains.add(headerInfo.getHeaderStrings()[i]);
        }
        ArrayList<String> sortedStrains = new ArrayList<String>(
                unsortedStrains);
        sortedStrains.retainAll(strainsToParse);
        Collections.sort(sortedStrains);
        int[] sortToUnsortIndexMapping = new int[sortedStrains.size()];
        for(int i = 0; i < sortToUnsortIndexMapping.length; i++)
        {
            int mappedIndex = unsortedStrains.indexOf(sortedStrains.get(i));
            sortToUnsortIndexMapping[i] = mappedIndex;
        }
        
        // make sure there are no duplicate strains
        if((new HashSet<String>(unsortedStrains)).size() != unsortedStrains.size())
        {
            throw new IOException(
                    "Found duplicate strain names in file header which is " +
                    "not allowed");
        }
        
        // OK now parse through the rest of the file and build the chromosome set
        int rowCount = 0;
        Map<Integer, List<MultiPartitionedInterval>> chromosomeStateMap =
            new TreeMap<Integer, List<MultiPartitionedInterval>>();
        String[] currRow;
        
        // variables needed to represent interval start
        short[] intervalStrainGroups = null;
        int intervalChromosomeNum = -1;
        long intervalStartBpPosition = -1;
        long intervalEndBpPosition = -1;
        
        while((currRow = this.characterDelimitedParser.parseCharacterDelimitedLine(bufferedReader)) != null)
        {
            int currChromosomeNum = this.getChromosomeNumber(
                    headerInfo,
                    currRow);
            long currBasePairPosition = this.getBasePairPosition(
                    headerInfo,
                    currRow);
            
            short[] currStrainGroups = new short[sortToUnsortIndexMapping.length];
            for(int i = 0; i < currStrainGroups.length; i++)
            {
                int unsortedIndex = sortToUnsortIndexMapping[i];
                int headerIndex = unsortedIndex + this.numAnnotationColumns;
                String strainGroupString = currRow[headerIndex];
                if(strainGroupString.equals("-"))
                {
                    // "-" means that we should skip the row
                    currStrainGroups = null;
                    break;
                }
                else
                {
                    currStrainGroups[i] = Short.parseShort(strainGroupString);
                }
            }
            
            if(currStrainGroups != null)
            {
                if(intervalStrainGroups == null)
                {
                    // we are at the start of a new interval
                    intervalStrainGroups = currStrainGroups;
                    intervalChromosomeNum = currChromosomeNum;
                    intervalStartBpPosition = currBasePairPosition;
                    intervalEndBpPosition = currBasePairPosition;
                }
                else if(this.flattenMatchingHaplotypes &&
                        Arrays.equals(intervalStrainGroups, currStrainGroups) &&
                        intervalChromosomeNum == currChromosomeNum)
                {
                    // since the curr values are consistent lets move the
                    // interval forward
                    intervalEndBpPosition = currBasePairPosition;
                }
                else
                {
                    // we've reached the end of an interval. close it out
                    this.addIntervalToMap(
                            chromosomeStateMap,
                            intervalStrainGroups,
                            intervalChromosomeNum,
                            intervalStartBpPosition,
                            intervalEndBpPosition);
                    
                    // start a new interval
                    intervalStrainGroups = currStrainGroups;
                    intervalChromosomeNum = currChromosomeNum;
                    intervalStartBpPosition = currBasePairPosition;
                    intervalEndBpPosition = currBasePairPosition;
                }
            }
            
            if(rowCount % 10000 == 0)
            {
                if(LOG.isLoggable(Level.FINE))
                {
                    LOG.fine("completed parsing row #: " + rowCount);
                }
                System.out.print('*');
            }
            rowCount++;
        }
        
        LOG.fine("finished parsing");
        
        // close out the final interval
        if(intervalStrainGroups != null)
        {
            this.addIntervalToMap(
                    chromosomeStateMap,
                    intervalStrainGroups,
                    intervalChromosomeNum,
                    intervalStartBpPosition,
                    intervalEndBpPosition);
        }
        
        // shrink the lists
        for(List<MultiPartitionedInterval> values :
            chromosomeStateMap.values())
        {
            ((ArrayList<MultiPartitionedInterval>)values).trimToSize();
        }
        
        return chromosomeStateMap;
    }
    
    private void addIntervalToMap(
            Map<Integer, List<MultiPartitionedInterval>> chromosomeStateMap,
            short[] intervalStrainGroups,
            int intervalChromosomeNum,
            long intervalStartBpPosition,
            long intervalEndBpPosition)
    {
        List<MultiPartitionedInterval> intervalList =
            chromosomeStateMap.get(intervalChromosomeNum);
        if(intervalList == null)
        {
            intervalList = new ArrayList<MultiPartitionedInterval>();
            chromosomeStateMap.put(intervalChromosomeNum, intervalList);
        }
        
        intervalList.add(new MultiPartitionedInterval(
                intervalChromosomeNum,
                intervalStartBpPosition,
                1 + intervalEndBpPosition - intervalStartBpPosition,
                intervalStrainGroups));
    }

    /**
     * Parse the 1st chromosome
     * @param bufferedReader
     *          the reader that we're parsing
     * @return
     *          the map of chromosomes with haplotypes
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public int parseFirstChromosome(BufferedReader bufferedReader) throws IOException
    {
        HeaderInfo headerInfo = this.parseHeaderInfoFromReader(bufferedReader);
        
        String[] firstRow =
            this.characterDelimitedParser.parseCharacterDelimitedLine(
                    bufferedReader);
        int chromosomeNum = this.getChromosomeNumber(
                headerInfo,
                firstRow);
        
        return chromosomeNum;
    }
}
