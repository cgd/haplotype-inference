package org.jax.haplotype.phylogeny.inference;

import org.jax.geneticutil.data.SimpleBasePairInterval;

/**
 * A {@link SimpleBasePairInterval} that also has a {@link SdpInclusionHierarchy}
 * @author Keith Sheppard
 */
/*package-protected*/ class SdpInclusionInterval extends SimpleBasePairInterval
{
    /**
     * every {@link java.io.Serializable} should have one of these
     */
    private static final long serialVersionUID = 3109196877163974973L;
    
    private final SdpInclusionHierarchy sdpInclusionHierarchy;

    /**
     * Constructor
     * @param chromosomeNumber
     *          see {@link #getChromosomeNumber()}
     * @param startInBasePairs
     *          see {@link #getStartInBasePairs()}
     * @param extentInBasePairs
     *          see {@link #getExtentInBasePairs()}
     * @param sdpInclusionHierarchy
     *          see {@link #getSdpInclusionHierarchy()}
     */
    public SdpInclusionInterval(
            int chromosomeNumber,
            long startInBasePairs,
            long extentInBasePairs,
            SdpInclusionHierarchy sdpInclusionHierarchy)
    {
        super(chromosomeNumber, startInBasePairs, extentInBasePairs);
        this.sdpInclusionHierarchy = sdpInclusionHierarchy;
    }
    
    /**
     * Getter for the SDP inclusion hierarchy
     * @return the sdpInclusionHierarchy
     */
    public SdpInclusionHierarchy getSdpInclusionHierarchy()
    {
        return this.sdpInclusionHierarchy;
    }
}
