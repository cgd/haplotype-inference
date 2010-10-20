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

package org.jax.haplotype;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFileChooser;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.haplotype.data.ChromosomeDataSource;
import org.jax.haplotype.data.CommaSeparatedChromosomeDataSource;
import org.jax.haplotype.inference.HaplotypeEstimator;
import org.jax.haplotype.inference.IntervalScanningHaplotypeEstimator;

/**
 * A main class for doing some haplotype inference
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HAM implements Runnable
{
    private static final Logger LOG =
        Logger.getLogger(HAM.class.getName());
    
    private final HaplotypeEstimator haplotypeEstimator =
        new IntervalScanningHaplotypeEstimator(3, 5);
    
    /**
     * Constructor
     */
    public HAM()
    {
        
    }
    
    /**
     * {@inheritDoc}
     */
    public void run()
    {
        JFileChooser fileChooser = new JFileChooser();
        int userSelection = fileChooser.showOpenDialog(null);
        if(userSelection == JFileChooser.APPROVE_OPTION)
        {
            File selectedFile = fileChooser.getSelectedFile();
            ChromosomeDataSource chromoDataSource =
                new CommaSeparatedChromosomeDataSource(
                        selectedFile,
                        -1);
            try
            {
                String[] sortedStrains = chromoDataSource.getAvailableStrains().toArray(
                        new String[0]);
                Arrays.sort(sortedStrains);
                
                for(String currStrainName: sortedStrains)
                {
                    System.out.println("Parsed Chromosome:");
                    System.out.println("  Strain: " + currStrainName);
                }
                
                List<PartitionedInterval> haplotypes =
                    this.haplotypeEstimator.estimateHaplotypeBlocks(
                            chromoDataSource.getSdpInputStream(sortedStrains),
                            chromoDataSource.getSnpPositionInputStream());
                for(PartitionedInterval currHaplotype: haplotypes)
                {
                    System.out.println("Estimated Haplotype: ");
                    System.out.println("  haplotype chromo count: " +
                            currHaplotype.getStrainBitSet().cardinality());
                    System.out.println("  haplotype start SNP:    " +
                            currHaplotype.getStartInBasePairs());
                    System.out.println("  haplotype extent:       " +
                            currHaplotype.getExtentInBasePairs());
                }
            }
            catch(Exception ex)
            {
                LOG.log(Level.SEVERE,
                        "failed to estimate haplotypes",
                        ex);
            }
        }
        else
        {
            System.out.println("user doesn't want to open the file");
        }
    }

    /**
     * The main HAM application entry point.
     * @param args
     *          don't use these
     */
    public static void main(String[] args)
    {
        Thread thread = new Thread(new HAM());
        thread.start();
    }

}
