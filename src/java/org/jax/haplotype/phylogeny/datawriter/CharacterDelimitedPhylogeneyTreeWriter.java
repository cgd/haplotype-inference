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

package org.jax.haplotype.phylogeny.datawriter;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.haplotype.phylogeny.data.PhylogenyTestResult;
import org.jax.haplotype.phylogeny.data.PhylogenyTreeNode;
import org.jax.util.io.CharacterDelimitedParser;

/**
 * Write out the trees in character delimited format
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CharacterDelimitedPhylogeneyTreeWriter
{
    private final CharacterDelimitedParser parser;
    
    private static final double EPSILON = 0.00001;
    
    /**
     * Constructor
     * @param parser
     *          the parser to use
     */
    public CharacterDelimitedPhylogeneyTreeWriter(
            CharacterDelimitedParser parser)
    {
        this.parser = parser;
    }
    
    /**
     * Same as calling
     * {@link #writePhylogenyIntervalsToStream(OutputStream, Collection, String, boolean)}
     * with {@code simlifyTrees} set to {@code false}
     * @param outputStream
     *          the output stream
     * @param phylogenies
     *          the phylogeny trees
     * @param headerComment
     *          the header comment. this comment does not have to have a
     *          '#' prefix since this function will add one
     * @throws IOException
     *          if we catch an {@link java.io.IOException} from the
     *          given output stream
     */
    public void writePhylogenyIntervalsToStream(
            OutputStream outputStream,
            Collection<PhylogenyInterval> phylogenies,
            String headerComment)
    throws IOException
    {
        this.writePhylogenyIntervalsToStream(
                outputStream,
                phylogenies,
                headerComment,
                false);
    }
    
    /**
     * Write the given phylogenies to the given output stream and tack
     * on the given header comment
     * @param outputStream
     *          the output stream
     * @param phylogenies
     *          the phylogeny trees
     * @param headerComment
     *          the header comment. this comment does not have to have a
     *          '#' prefix since this function will add one
     * @param simplifyTrees
     *          if true we should simplify trees using
     *          {@link PhylogenyTreeNode#removeNonBranchingInteriorNodes()} and
     *          {@link PhylogenyTreeNode#resolveToSingleStrainLeafNodes(double)}
     *          before writing the tree to file
     * @throws IOException
     *          if we catch an {@link java.io.IOException} from the
     *          given output stream
     */
    public void writePhylogenyIntervalsToStream(
            OutputStream outputStream,
            Collection<PhylogenyInterval> phylogenies,
            String headerComment,
            boolean simplifyTrees)
    throws IOException
    {
        // make sure to buffer the output
        if(!(outputStream instanceof BufferedOutputStream))
        {
            outputStream = new BufferedOutputStream(outputStream);
        }
        PrintStream printStream = new PrintStream(
                outputStream);
        
        // write out the header comment
        if(headerComment != null && headerComment.length() > 0)
        {
            this.parser.writeHeaderComment(
                    headerComment,
                    printStream);
        }
        
        this.parser.writeCharacterDelimitedRow(new String[] {
                "chromosomeNumber",
                "intervalStartInBasePairs",
                "intervalExtentInBasePairs",
                "newickPerfectPhylogeny"},
                printStream);
        
        for(PhylogenyInterval currPhylogeny: phylogenies)
        {
            PhylogenyTreeNode tree = currPhylogeny.getPhylogeny();
            
            if(simplifyTrees)
            {
                // simplify the tree so that analysis is easier
                tree = tree.resolveToSingleStrainLeafNodes(EPSILON);
                tree = tree.removeNonBranchingInteriorNodes();
            }
            
            BasePairInterval interval =
                currPhylogeny.getInterval();
            
            // print the tree
            String[] row = new String[] {
                    Integer.toString(interval.getChromosomeNumber()),
                    Long.toString(interval.getStartInBasePairs()),
                    Long.toString(interval.getExtentInBasePairs()),
                    tree.toNewickFormat()};
            this.parser.writeCharacterDelimitedRow(row, printStream);
        }
        
        printStream.flush();
    }
    
    /**
     * Write the phylogeny test results to stream
     * @param outputStream
     *          the output stream
     * @param testResults
     *          the test results
     */
    public void writePhylogenyTestResultsToStream(
            OutputStream outputStream,
            Collection<PhylogenyTestResult> testResults)
    {
        this.writePhylogenyTestResultsToStream(
                outputStream,
                testResults,
                null);
    }
    
    /**
     * Write the phylogeny test results to stream
     * @param outputStream
     *          the output stream
     * @param testResults
     *          the test results
     * @param headerComment
     *          the header comment
     */
    public void writePhylogenyTestResultsToStream(
            OutputStream outputStream,
            Collection<PhylogenyTestResult> testResults,
            String headerComment)
    {
        // make sure to buffer the output
        if(!(outputStream instanceof BufferedOutputStream))
        {
            outputStream = new BufferedOutputStream(outputStream);
        }
        PrintStream printStream = new PrintStream(
                outputStream);
        
        // write out the header comment
        if(headerComment != null && headerComment.length() > 0)
        {
            this.parser.writeHeaderComment(
                    headerComment,
                    printStream);
        }
        
        this.parser.writeCharacterDelimitedRow(new String[] {
                "chromosomeNumber",
                "intervalStartInBasePairs",
                "intervalExtentInBasePairs",
                "newickPerfectPhylogeny",
                "treePValue"},
                printStream);
        
        for(PhylogenyTestResult testResult: testResults)
        {
            BasePairInterval interval =
                testResult.getPhylogenyInterval().getInterval();
            PhylogenyTreeNode tree =
                testResult.getPhylogenyInterval().getPhylogeny();
            
            String[] row = new String[] {
                    Integer.toString(interval.getChromosomeNumber()),
                    Long.toString(interval.getStartInBasePairs()),
                    Long.toString(interval.getExtentInBasePairs()),
                    tree.toNewickFormat(),
                    Double.toString(testResult.getPValue())};
            this.parser.writeCharacterDelimitedRow(row, printStream);
        }
    }
}
