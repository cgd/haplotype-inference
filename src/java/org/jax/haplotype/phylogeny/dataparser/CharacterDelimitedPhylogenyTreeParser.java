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

package org.jax.haplotype.phylogeny.dataparser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.SimpleBasePairInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.io.CharacterDelimitedParser;
import org.jax.util.io.IllegalFormatException;

/**
 * A parser for phylogeny data.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CharacterDelimitedPhylogenyTreeParser
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            CharacterDelimitedPhylogenyTreeParser.class.getName());
    
    private final CharacterDelimitedParser parser;

    /**
     * Constructor
     * @param parser
     *          the parser to use
     */
    public CharacterDelimitedPhylogenyTreeParser(CharacterDelimitedParser parser)
    {
        this.parser = parser;
    }
    
    /**
     * Parse the data from the given formatted input stream and generate
     * phylogeny intervals as output
     * @param is
     *          the input stream to parse through
     * @return
     *          all of the intervals mapped by chromosome number
     * @throws IOException
     */
    public Map<Integer, List<PhylogenyInterval>> parsePhylogenyIntervalsFromStream(
            InputStream is)
            throws IOException
    {
        BufferedReader bufferedReader = new BufferedReader(
                new InputStreamReader(is));
        
        String currLine = bufferedReader.readLine();
        while(currLine != null && this.parser.isCommentOrEmpty(currLine))
        {
            currLine = bufferedReader.readLine();
        }
        
        if(currLine == null)
        {
            LOG.warning("did not find any phylogeny intervals in stream");
            return Collections.emptyMap();
        }
        else
        {
            Map<Integer, List<PhylogenyInterval>> intervalMap =
                new HashMap<Integer, List<PhylogenyInterval>>();
            
            String[] currLineTokens = this.parser.parseCharacterDelimitedLine(
                    currLine);
            do
            {
                int chromosomeNumber = Integer.parseInt(currLineTokens[0]);
                long startBp = Long.parseLong(currLineTokens[1]);
                long extentBp = Long.parseLong(currLineTokens[2]);
                String newickTreeString = currLineTokens[3];
                
                try
                {
                    PhylogenyTreeNode tree = PhylogenyTreeNode.fromNewickFormat(
                            newickTreeString);
                    BasePairInterval snpInterval =
                        new SimpleBasePairInterval(
                                chromosomeNumber,
                                startBp,
                                extentBp);
                    PhylogenyInterval treeInterval = new PhylogenyInterval(
                            tree,
                            snpInterval);
                    
                    List<PhylogenyInterval> treeIntervalList = intervalMap.get(
                            chromosomeNumber);
                    if(treeIntervalList == null)
                    {
                        treeIntervalList = new ArrayList<PhylogenyInterval>();
                        intervalMap.put(chromosomeNumber, treeIntervalList);
                    }
                    treeIntervalList.add(treeInterval);
                }
                catch(IllegalFormatException ex)
                {
                    // TODO for JDK 6 throw this as a cause
                    LOG.log(Level.SEVERE,
                            "Bad format",
                            ex);
                    throw new IOException(ex.getMessage());
                }
            } while((currLineTokens = this.parser.parseCharacterDelimitedLine(bufferedReader)) != null);
            
            return intervalMap;
        }
    }
}
