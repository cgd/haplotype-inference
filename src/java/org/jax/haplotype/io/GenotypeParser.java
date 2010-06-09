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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.SnpType;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.util.io.CharacterDelimitedParser;

/**
 * Reads genotype information from a CSV file.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class GenotypeParser extends GenomicFlatFileParser
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            GenotypeParser.class.getName());
    
    /**
     * Constructor
     */
    public GenotypeParser()
    {
        super();
    }
    
    /**
     * Constructor
     * @param characterDelimitedParser
     *          the parser
     * @param numAnnotationColumns
     *          the number of annotation columns
     */
    public GenotypeParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns)
    {
        super(characterDelimitedParser, numAnnotationColumns);
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
    public GenotypeParser(
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition)
    {
        super(new CharacterDelimitedParser(),
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
    public GenotypeParser(
            CharacterDelimitedParser characterDelimitedParser,
            int numAnnotationColumns,
            int indexOfChromosomeId,
            int indexOfBasePairPosition)
    {
        super(characterDelimitedParser,
              numAnnotationColumns,
              indexOfChromosomeId,
              indexOfBasePairPosition);
    }
    
    /**
     * Parse the input stream into a set of chromosomes.
     * @param is
     *          the input stream that we're parsing
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromStream(InputStream is) throws IOException
    {
        return this.parseGenotypeFromStream(is, true, null);
    }
    
    /**
     * Parse the input stream into a set of chromosomes.
     * @param is
     *          the input stream that we're parsing
     * @param allowImputedSnps
     *          if true we'll read in the inputed SNPs, otherwise we
     *          discard them
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromStream(
            InputStream is,
            boolean allowImputedSnps) throws IOException
    {
        return this.parseGenotypeFromStream(is, allowImputedSnps, null);
    }
    
    /**
     * Parse the input stream into a set of chromosomes.
     * @param inputStream
     *          the input stream that we're parsing
     * @param allowImputedSnps
     *          if true we'll read in the inputed SNPs, otherwise we
     *          discard them
     * @param strainsToParse
     *          the names of the strains we should parse... if null parse
     *          everything
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromStream(
            InputStream inputStream,
            boolean allowImputedSnps,
            Set<String> strainsToParse) throws IOException
    {
        return this.parseGenotypeFromReader(
                new BufferedReader(new InputStreamReader(inputStream)),
                allowImputedSnps,
                strainsToParse);
    }
    
    /**
     * Parse the reader into a set of chromosomes.
     * @param reader
     *          the reader that we're parsing
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromReader(BufferedReader reader) throws IOException
    {
        return this.parseGenotypeFromReader(reader, true, null);
    }
    
    /**
     * Parse the reader into a set of chromosomes.
     * @param reader
     *          the reader that we're parsing
     * @param allowImputedSnps
     *          if true we'll read in the inputed SNPs, otherwise we
     *          discard them
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromReader(
            BufferedReader reader,
            boolean allowImputedSnps) throws IOException
    {
        return this.parseGenotypeFromReader(reader, allowImputedSnps, null);
    }
    
    /**
     * Parse the reader into a set of chromosomes.
     * @param bufferedReader
     *          the reader that we're parsing
     * @param allowImputedSnps
     *          if true we'll read in the inputed SNPs, otherwise we
     *          discard them
     * @param strainsToParse
     *          the names of the strains we should parse... if null parse
     *          everything
     * @return
     *          the set of chromosomes that we parsed through
     * @throws IOException
     *          if we get any exceptions from the input stream
     */
    public Set<StrainChromosome> parseGenotypeFromReader(
            BufferedReader bufferedReader,
            boolean allowImputedSnps,
            Set<String> strainsToParse) throws IOException
    {
        System.out.print("parsing genotype data: ");
        
        HeaderInfo headerInfo = this.parseHeaderInfoFromReader(bufferedReader);
        
        // OK now parse through the rest of the file and build the chromosome set
        int rowCount = 0;
        Map<StrainChromosome, List<SingleNucleotidePolymorphism>> snpMap =
            new TreeMap<StrainChromosome, List<SingleNucleotidePolymorphism>>();
        String[] currRow;
        while((currRow = this.characterDelimitedParser.parseCharacterDelimitedLine(bufferedReader)) != null)
        {
            int chromosomeNum = this.getChromosomeNumber(
                    headerInfo,
                    currRow);
            long currBasePairPosition = this.getBasePairPosition(
                    headerInfo,
                    currRow);
            
            for(int i = this.numAnnotationColumns;
                i < headerInfo.getRowLengthWithoutConflictColumn();
                i++)
            {
                String currStrain = headerInfo.getHeaderStrings()[i];
                if(strainsToParse == null || strainsToParse.contains(currStrain))
                {
                    // lower case means imputed
                    if(allowImputedSnps || Character.isUpperCase(currRow[i].charAt(0)))
                    {
                        // we need to get our hands on the SNP list for the right chromosome
                        StrainChromosome currChromosome = new StrainChromosome(
                                currStrain,
                                chromosomeNum);
                        List<SingleNucleotidePolymorphism> currChromosomeSNPList =
                            snpMap.get(currChromosome);
                        if(currChromosomeSNPList == null)
                        {
                            // couldn't find the SNP list for the chromosome so we'll make a new one
                            currChromosomeSNPList = new ArrayList<SingleNucleotidePolymorphism>();
                            snpMap.put(currChromosome, currChromosomeSNPList);
                        }
                        
                        // TODO get a real confidence value
                        // now we need to create the SNP
                        SingleNucleotidePolymorphism currSNP = new SingleNucleotidePolymorphism(
                                SnpType.snpStringToSNPEnum(currRow[i]),
                                currBasePairPosition);
                        currChromosomeSNPList.add(currSNP);
                    }
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
        
        // now we can turn our map into a set
        LOG.fine("finished parsing converting the result to arrays");
        Set<StrainChromosome> chromosomeSet =
            new HashSet<StrainChromosome>();
        Iterator<Map.Entry<StrainChromosome, List<SingleNucleotidePolymorphism>>> iter = snpMap.entrySet().iterator();
        while(iter.hasNext())
        {
            // Get the next element and remove it right away so it can be garbage collected
            Map.Entry<StrainChromosome, List<SingleNucleotidePolymorphism>> currChromosomeEntry = iter.next();
            iter.remove();
            
            StrainChromosome currChromosome = currChromosomeEntry.getKey();
            List<SingleNucleotidePolymorphism> currSNPList = currChromosomeEntry.getValue();
            currChromosome.setSingleNucleotidePolymorphisms(
                    currSNPList.toArray(new SingleNucleotidePolymorphism[currSNPList.size()]));
            chromosomeSet.add(currChromosome);
        }
        System.out.println();
        
        return chromosomeSet;
    }
}
