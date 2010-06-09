package org.jax.haplotype.io;

import java.io.IOException;
import java.util.BitSet;
import java.util.LinkedList;
import java.util.List;

/**
 * An SDP input stream that lets you do stuff similar in concept to a
 * {@link java.io.FilterInputStream}
 * @author Keith Sheppard
 */
public class BufferedSdpInputStream implements SdpInputStream
{
    private final SdpInputStream delegateSdpInputStream;
    
    private List<BitSet> activeBuffer = new LinkedList<BitSet>();
    
    private List<BitSet> markBuffer = new LinkedList<BitSet>();
    
    /**
     * Constructor
     * @param delegateSdpInputStream
     *          
     */
    public BufferedSdpInputStream(SdpInputStream delegateSdpInputStream)
    {
        this.delegateSdpInputStream = delegateSdpInputStream;
    }
    
    /**
     * Similar in concept to the {@link java.io.FilterInputStream#mark(int)}
     * function. Marks a point in the stream that we can later
     * {@link #reset()} to.
     */
    public void mark()
    {
        this.markBuffer = new LinkedList<BitSet>();
    }
    
    /**
     * Reset the stream to the last mark position
     */
    public void reset()
    {
        if(this.markBuffer == null)
        {
            throw new IllegalStateException(
                    "no current mark set");
        }
        
        if(this.activeBuffer != null)
        {
            this.markBuffer.addAll(this.activeBuffer);
        }
        
        this.activeBuffer = this.markBuffer;
        this.markBuffer = null;
    }
    
    /**
     * Push the given SDP bit set into the next position in the stream
     * @param sdp
     *          the SDP to push
     */
    public void push(BitSet sdp)
    {
        if(this.activeBuffer == null)
        {
            this.activeBuffer = new LinkedList<BitSet>();
        }
        
        this.activeBuffer.add(0, sdp);
    }

    /**
     * {@inheritDoc}
     */
    public BitSet getNextSdp() throws IOException
    {
        final BitSet nextSdp;
        if(this.activeBuffer != null)
        {
            nextSdp = this.activeBuffer.remove(0);
            if(this.activeBuffer.isEmpty())
            {
                this.activeBuffer = null;
            }
        }
        else
        {
            nextSdp = this.delegateSdpInputStream.getNextSdp();
        }
        
        if(this.markBuffer != null)
        {
            this.markBuffer.add(nextSdp);
        }
        
        return nextSdp;
    }

    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return this.delegateSdpInputStream.getReadDirection();
    }

    /**
     * {@inheritDoc}
     */
    public long getSdpCount() throws IOException
    {
        return this.delegateSdpInputStream.getSdpCount();
    }

    /**
     * {@inheritDoc}
     */
    public String[] getSdpStrainNames() throws IOException
    {
        return this.delegateSdpInputStream.getSdpStrainNames();
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSdp() throws IOException
    {
        return this.delegateSdpInputStream.hasNextSdp() ||
               !this.activeBuffer.isEmpty();
    }
}
