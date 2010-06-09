package org.jax.haplotype.phylogeny.inference;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Comparator;
import java.util.List;

import org.jax.haplotype.phylogeny.data.PhylogenyInterval;
import org.jax.util.datastructure.SetUtilities;

/**
 * A recursive structure representing the SDP inclusion hierarchy for a
 * perfect phylogeny. This is sort of a different way of looking at the
 * similar kind of data to {@link PhylogenyInterval}.
 * @author Keith Sheppard
 */
/*package-protected*/ class SdpInclusionHierarchy
{
    /**
     * A comparator that only looks at the SDP bits. It's capable of taking
     * an {@link SdpInclusionHierarchy} and/or a {@link BitSet} as input
     * @see SetUtilities#BIT_SET_COMPARATOR
     */
    public static final Comparator<Object> SDP_ONLY_COMPARATOR =
        new Comparator<Object>()
        {
            /**
             * {@inheritDoc}
             */
            public int compare(Object object1, Object object2)
            {
                final BitSet sdpBits1;
                if(object1 instanceof BitSet)
                {
                    sdpBits1 = (BitSet)object1;
                }
                else if(object1 instanceof SdpInclusionHierarchy)
                {
                    SdpInclusionHierarchy inclusionHierarchy1 =
                        (SdpInclusionHierarchy)object1;
                    sdpBits1 = inclusionHierarchy1.getSdpBits();
                }
                else
                {
                    throw new IllegalArgumentException(
                            "argument 1 must be a BitSet or an " +
                            "SdpInclusionHierarchy");
                }
                
                final BitSet sdpBits2;
                if(object2 instanceof BitSet)
                {
                    sdpBits2 = (BitSet)object2;
                }
                else if(object2 instanceof SdpInclusionHierarchy)
                {
                    SdpInclusionHierarchy inclusionHierarchy2 =
                        (SdpInclusionHierarchy)object2;
                    sdpBits2 = inclusionHierarchy2.getSdpBits();
                }
                else
                {
                    throw new IllegalArgumentException(
                            "argument 2 must be a BitSet or an " +
                            "SdpInclusionHierarchy");
                }
                
                return SetUtilities.BIT_SET_COMPARATOR.compare(
                        sdpBits1,
                        sdpBits2);
            }
        };
    
    private final BitSet sdpBits;
    
    private final List<SdpInclusionHierarchy> children;

    /**
     * Constructor
     * @param sdpBits
     *          {@link #getSdpBits()}
     */
    public SdpInclusionHierarchy(BitSet sdpBits)
    {
        this(sdpBits, new ArrayList<SdpInclusionHierarchy>());
    }
    
    /**
     * Constructor
     * @param sdpBits
     *          see {@link #getSdpBits()}
     * @param children
     *          see {@link #getChildren()}
     */
    public SdpInclusionHierarchy(
            BitSet sdpBits,
            List<SdpInclusionHierarchy> children)
    {
        this.sdpBits = sdpBits;
        this.children = children;
    }
    
    /**
     * Getter for the SDP bits
     * @return the sdpBits
     */
    public BitSet getSdpBits()
    {
        return this.sdpBits;
    }
    
    /**
     * Getter for the children
     * @return the children
     */
    public List<SdpInclusionHierarchy> getChildren()
    {
        return this.children;
    }
}
