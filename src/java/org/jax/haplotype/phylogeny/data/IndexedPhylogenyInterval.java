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

import java.io.Serializable;

import org.jax.geneticutil.data.IndexedSnpInterval;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class IndexedPhylogenyInterval implements Serializable
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = 6181753762556468516L;

    private PhylogenyTreeNode phylogeny;
    
    private IndexedSnpInterval indexedInterval;
    
    /**
     * Constructor
     * @param indexedInterval
     *          see {@link #getIndexedInterval()}
     * @param phylogeny
     *          see {@link #getPhylogeny()}
     */
    public IndexedPhylogenyInterval(
            PhylogenyTreeNode phylogeny,
            IndexedSnpInterval indexedInterval)
    {
        this.indexedInterval = indexedInterval;
        this.phylogeny = phylogeny;
    }
    
    /**
     * Getter for the phylogeny
     * @return the phylogeny
     */
    public PhylogenyTreeNode getPhylogeny()
    {
        return this.phylogeny;
    }
    
    /**
     * Getter for the indexed interval
     * @return the indexed interval
     */
    public IndexedSnpInterval getIndexedInterval()
    {
        return this.indexedInterval;
    }
}
