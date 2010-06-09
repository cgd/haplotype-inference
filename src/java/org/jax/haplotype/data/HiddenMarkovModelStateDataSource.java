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
