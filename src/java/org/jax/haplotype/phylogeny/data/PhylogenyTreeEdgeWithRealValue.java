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

package org.jax.haplotype.phylogeny.data;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import org.jax.util.math.Function;

/**
 * A tree edge with a real value associated with it
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyTreeEdgeWithRealValue extends PhylogenyTreeEdge
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 5628687239159596156L;
    private double realValue;
    
    /**
     * Constructor
     * @param sdpBits
     *          see {@link #getSdpBits()}
     * @param node
     *          see {@link #getNode()}
     * @param edgeLength
     *          see {@link #getEdgeLength()}
     * @param realValue
     *          see {@link #getRealValue()}
     */
    public PhylogenyTreeEdgeWithRealValue(
            BitSet sdpBits,
            PhylogenyTreeNode node,
            double edgeLength,
            double realValue)
    {
        super(sdpBits, node, edgeLength);
        this.realValue = realValue;
    }

    /**
     * Getter for the value
     * @return the realValue
     */
    public double getRealValue()
    {
        return this.realValue;
    }
    
    /**
     * Setter for the value
     * @param realValue the realValue to set
     */
    public void setRealValue(double realValue)
    {
        this.realValue = realValue;
    }
    
    /**
     * Use a function to transform the real value in a tree
     * @param node
     *          the node to transform
     * @param function
     *          the function
     * @return
     *          the node with all {@link PhylogenyTreeEdgeWithRealValue} edges
     *          transformed by the given function
     */
    public static PhylogenyTreeNode transform(
            PhylogenyTreeNode node,
            Function<Double, ? extends Number> function)
    {
        List<PhylogenyTreeEdge> childEdges =
            new ArrayList<PhylogenyTreeEdge>(node.getChildEdges().size());
        for(PhylogenyTreeEdge edge: node.getChildEdges())
        {
            PhylogenyTreeNode transNode = transform(edge.getNode(), function);
            if(edge instanceof PhylogenyTreeEdgeWithRealValue)
            {
                PhylogenyTreeEdgeWithRealValue edgeWithValue =
                    (PhylogenyTreeEdgeWithRealValue)edge;
                Number newValue = function.evaluate(
                        Double.valueOf(edgeWithValue.getRealValue()));
                
                PhylogenyTreeEdgeWithRealValue newEdge =
                    new PhylogenyTreeEdgeWithRealValue(
                            edgeWithValue.getSdpBits(),
                            transNode,
                            edgeWithValue.getEdgeLength(),
                            newValue.doubleValue());
                childEdges.add(newEdge);
            }
        }
        
        PhylogenyTreeNode newNode = new PhylogenyTreeNode(
                childEdges,
                node.getStrains());
        return newNode;
    }
    
    /**
     * Get the smallest of all of the edges
     * @param node
     *          the node to extract the minimum edge from
     * @return
     *          the minimum edge
     */
    public static PhylogenyTreeEdgeWithRealValue getEdgeWithMininumValue(PhylogenyTreeNode node)
    {
        return getEdgeWithMinimumValueRecursive(
                node,
                null);
    }
    
    private static PhylogenyTreeEdgeWithRealValue getEdgeWithMinimumValueRecursive(
            PhylogenyTreeNode node,
            PhylogenyTreeEdgeWithRealValue minEdge)
    {
        for(PhylogenyTreeEdge currEdge: node.getChildEdges())
        {
            if(currEdge instanceof PhylogenyTreeEdgeWithRealValue)
            {
                PhylogenyTreeEdgeWithRealValue currEdgeWithReal =
                    (PhylogenyTreeEdgeWithRealValue)currEdge;
                if(minEdge == null ||
                   currEdgeWithReal.getRealValue() < minEdge.getRealValue())
                {
                    minEdge = currEdgeWithReal;
                }
            }
        }
        
        return minEdge;
    }
    
    /**
     * Get the biggest of all of the edges
     * @param node
     *          the node to extract the minimum edge from
     * @return
     *          the minimum edge
     */
    public static PhylogenyTreeEdgeWithRealValue getEdgeWithMaximumValue(PhylogenyTreeNode node)
    {
        return getEdgeWithMaximumValueRecursive(
                node,
                null);
    }
    
    private static PhylogenyTreeEdgeWithRealValue getEdgeWithMaximumValueRecursive(
            PhylogenyTreeNode node,
            PhylogenyTreeEdgeWithRealValue maxEdge)
    {
        for(PhylogenyTreeEdge currEdge: node.getChildEdges())
        {
            if(currEdge instanceof PhylogenyTreeEdgeWithRealValue)
            {
                PhylogenyTreeEdgeWithRealValue currEdgeWithReal =
                    (PhylogenyTreeEdgeWithRealValue)currEdge;
                if(maxEdge == null ||
                   currEdgeWithReal.getRealValue() > maxEdge.getRealValue())
                {
                    maxEdge = currEdgeWithReal;
                }
            }
        }
        
        return maxEdge;
    }
}
