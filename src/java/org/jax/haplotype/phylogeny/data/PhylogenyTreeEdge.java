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

package org.jax.haplotype.phylogeny.data;

import java.io.Serializable;
import java.util.BitSet;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.jax.util.ObjectUtil;

/**
 * An edge of a perfect phylogeny tree
 * @see PhylogenyTreeNode
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class PhylogenyTreeEdge implements Cloneable, Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -7237477674690794882L;

    private static final Logger LOG = Logger.getLogger(
            PhylogenyTreeEdge.class.getName());
    
    private BitSet sdpBits;
    
    private PhylogenyTreeNode node;
    
    private double edgeLength;
    
    /**
     * Constructor
     */
    public PhylogenyTreeEdge()
    {
        this(null, null, 1.0);
    }
    
    /**
     * Constructor
     * @param sdpBits
     *          see {@link #getSdpBits()}
     * @param node
     *          see {@link #getNode()}
     * @param edgeLength
     *          see {@link #getEdgeLength()}
     */
    public PhylogenyTreeEdge(
            BitSet sdpBits,
            PhylogenyTreeNode node,
            double edgeLength)
    {
        this.sdpBits = sdpBits;
        this.node = node;
        this.edgeLength = edgeLength;
    }
    
    /**
     * Getter for the edge length. No units are implied at this level (you
     * should track units somewhere else)
     * @return the edgeLength
     */
    public double getEdgeLength()
    {
        return this.edgeLength;
    }
    
    /**
     * Setter for the edge length
     * @see #getEdgeLength()
     * @param edgeLength the edgeLength to set
     */
    public void setEdgeLength(double edgeLength)
    {
        this.edgeLength = edgeLength;
    }

    /**
     * Getter for the Strain Distribution pattern for this edge represented as
     * bit flags
     * @return
     *          the sdpBits
     */
    public BitSet getSdpBits()
    {
        return this.sdpBits;
    }
    
    /**
     * Setter for the SDP bits
     * @see #getSdpBits()
     * @param sdpBits the sdpBits to set
     */
    public void setSdpBits(BitSet sdpBits)
    {
        this.sdpBits = sdpBits;
    }
    
    /**
     * Getter for the child node at the end of this edge
     * @return the node
     */
    public PhylogenyTreeNode getNode()
    {
        return this.node;
    }
    
    /**
     * Setter for the child node at the end of this edge
     * @param node the node to set
     */
    public void setNode(PhylogenyTreeNode node)
    {
        this.node = node;
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    protected PhylogenyTreeEdge clone()
    {
        try
        {
            PhylogenyTreeEdge newEdge = (PhylogenyTreeEdge)super.clone();
            
            newEdge.node = newEdge.node.clone();
            return newEdge;
        }
        catch(CloneNotSupportedException ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to clone edge. this should never happen",
                    ex);
            return null;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public boolean equals(Object obj)
    {
        if(obj instanceof PhylogenyTreeEdge)
        {
            PhylogenyTreeEdge otherEdge = (PhylogenyTreeEdge)obj;
            return this.edgeLength == otherEdge.edgeLength &&
//                   ObjectUtil.areEqual(this.sdpBits, otherEdge.sdpBits) &&
                   ObjectUtil.areEqual(this.node, otherEdge.node);
        }
        else
        {
            return false;
        }
    }
    
    /**
     * {@inheritDoc}
     */
    @Override
    public int hashCode()
    {
        return this.node.hashCode();
    }
}
