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

package org.jax.ham.dataparser;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

import javax.swing.JFileChooser;

import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.GenotypeParser;

/**
 * for testing the parser
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommaSeparatedGenotypeParserTestMain
{
    /**
     * the main test method
     * @param args
     *          don't care about this
     * @throws IOException
     *          if we have io problems
     */
    public static void main(String[] args) throws IOException
    {
        GenotypeParser parser = new GenotypeParser();
        JFileChooser fileChooser = new JFileChooser();
        fileChooser.setMultiSelectionEnabled(true);
        int userSelection = fileChooser.showOpenDialog(null);
        if(userSelection == JFileChooser.APPROVE_OPTION)
        {
            File[] selectedFiles = fileChooser.getSelectedFiles();
            for(File selectedFile: selectedFiles)
            {
                FileInputStream fis = new FileInputStream(selectedFile);
                Set<StrainChromosome> chromoSet = parser.parseGenotypeFromStream(fis);
                
                StrainChromosome firstChromo = chromoSet.iterator().next();
                SingleNucleotidePolymorphism[] firstSnps =
                    firstChromo.getSingleNucleotidePolymorphisms();
                Set<SingleNucleotidePolymorphism> snpSet =
                    new HashSet<SingleNucleotidePolymorphism>();
                for(int i = 0; i < firstSnps.length; i++)
                {
                    for(StrainChromosome currChromosome: chromoSet)
                    {
                        snpSet.add(
                                currChromosome.getSingleNucleotidePolymorphisms()[i]);
                    }
                    
                    if(snpSet.size() == 0)
                    {
                        throw new IOException("0 snps at index " + i);
                    }
                    else if(snpSet.size() == 1)
                    {
                        System.out.println("only one snp type at " + i);
                    }
                    else if(snpSet.size() > 2)
                    {
                        throw new IOException("more than two snp types at: " + i);
                    }
                    
                    snpSet.clear();
                }
//                for(StrainChromosome currChromosome: chromoSet)
//                {
//                    System.out.println("Parsed Chromosome:");
//                    System.out.println("  Strain: " + currChromosome.getStrainName());
//                    System.out.println("  Number: " + currChromosome.getChromosomeNumber());
//                    
//                    for(SingleNucleotidePolymorphism currSNP: currChromosome.getSingleNucleotidePolymorphisms())
//                    {
//                        System.out.println(currSNP.getSnpType());
//                    }
//                }
            }
        }
        else
        {
            System.out.println("user doesn't want to open the file");
        }
    }
}
