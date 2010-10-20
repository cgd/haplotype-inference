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

package org.jax.haplotype.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.jax.geneticutil.data.MultiPartitionedInterval;
import org.jax.haplotype.io.HiddenMarkovModelStateParser;
import org.jax.util.datastructure.SequenceUtilities;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class HiddenMarkovModelStateDataSource implements
        MultiGroupHaplotypeDataSource
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -5814023098008501759L;

    private final HiddenMarkovModelStateParser parser;
    
    private final Map<Integer, File> chromosomeToCsvFileMap =
        new HashMap<Integer, File>();
    
    private final String name;
    
    /**
     * Constructor
     * @param name
     *          the name of the data source
     * @param files
     *          the files
     * @param parser
     *          the parser to use
     */
    public HiddenMarkovModelStateDataSource(
            String name,
            Collection<File> files,
            HiddenMarkovModelStateParser parser)
    {
        this.name = name;
        this.parser = parser;
        this.initialize(files);
    }
    
    private void initialize(Collection<File> files)
    {
        for(File file: files)
        {
            try
            {
                BufferedReader reader = new BufferedReader(new FileReader(file));
                try
                {
                    int chromo = this.parser.parseFirstChromosome(reader);
                    
                    if(this.chromosomeToCsvFileMap.put(chromo, file) != null)
                    {
                        throw new IOException(
                                "Error: two of the input files have " +
                                "chromosome number: " + chromo);
                    }
                }
                finally
                {
                    reader.close();
                }
            }
            catch(IOException ex)
            {
                throw new RuntimeException(ex);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public int[] getAvailableChromosomes()
    {
        int[] chrs = SequenceUtilities.toIntArray(new ArrayList<Integer>(
                this.chromosomeToCsvFileMap.keySet()));
        Arrays.sort(chrs);
        return chrs;
    }

    /**
     * {@inheritDoc}
     */
    public Set<String> getAvailableStrains()
    {
        if(this.chromosomeToCsvFileMap.isEmpty())
        {
            return Collections.emptySet();
        }
        else
        {
            File anyFile = this.chromosomeToCsvFileMap.values().iterator().next();
            
            try
            {
                FileInputStream fis = new FileInputStream(anyFile);
                
                try
                {
                    Set<String> availableStrains =
                        this.parser.parseAvailableStrainsFromStream(
                                fis);
                    return availableStrains;
                }
                finally
                {
                    fis.close();
                }
            }
            catch(IOException ex)
            {
                throw new RuntimeException(ex);
            }
        }
    }

    /**
     * {@inheritDoc}
     */
    public List<MultiPartitionedInterval> getHaplotypeData(
            int chromosome,
            Set<String> strainsToAccept)
    {
        File chrFile = this.chromosomeToCsvFileMap.get(chromosome);
        
        try
        {
            BufferedReader reader = new BufferedReader(new FileReader(chrFile));
            Map<Integer, List<MultiPartitionedInterval>> hmmStates =
                this.parser.parseHMMStatesFromReader(
                        reader,
                        strainsToAccept);
            
            try
            {
                if(hmmStates.size() >= 2)
                {
                    throw new RuntimeException(
                            "Error found more than one chromosome in: " +
                            chrFile.getAbsolutePath());
                }
                else if(hmmStates.isEmpty())
                {
                    throw new RuntimeException(
                            "Did not successfully parse: " +
                            chrFile.getAbsolutePath());
                }
                else
                {
                    return hmmStates.values().iterator().next();
                }
            }
            finally
            {
                reader.close();
            }
        }
        catch(IOException ex)
        {
            throw new RuntimeException(ex);
        }
    }

    /**
     * {@inheritDoc}
     */
    public String getName()
    {
        return this.name;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public String toString()
    {
        return this.getName();
    }
}
