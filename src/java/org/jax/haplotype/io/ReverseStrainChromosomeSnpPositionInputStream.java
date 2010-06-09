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

package org.jax.haplotype.io;

import java.io.IOException;

import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;

/**
 * A SNP position input stream wrapped around a {@link StrainChromosome}
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ReverseStrainChromosomeSnpPositionInputStream implements
        SnpPositionInputStream
{
    private final SingleNucleotidePolymorphism[] snps;
    
    private final int chromosomeNumber;
    
    private int index;
    
    /**
     * Constructor
     * @param strainChromosome
     *          the backing strain chromosome for this stream
     */
    public ReverseStrainChromosomeSnpPositionInputStream(
            StrainChromosome strainChromosome)
    {
        this(strainChromosome.getChromosomeNumber(),
             strainChromosome.getSingleNucleotidePolymorphisms());
    }
    
    /**
     * Constructor
     * @param chromosomeNumber
     *          the chromosome number for all snps
     * @param snps
     *          the backing snps for this stream
     */
    public ReverseStrainChromosomeSnpPositionInputStream(
            int chromosomeNumber,
            SingleNucleotidePolymorphism[] snps)
    {
        this.chromosomeNumber = chromosomeNumber;
        this.snps = snps;
        this.index = snps.length - 1;
    }

    /**
     * {@inheritDoc}
     */
    public long getNextSnpPositionInBasePairs() throws IOException
    {
        long nextSnpPosition = this.snps[this.index].getPositionInBasePairs();
        this.index--;
        return nextSnpPosition;
    }

    /**
     * {@inheritDoc}
     */
    public long getSnpCount() throws IOException
    {
        return this.snps.length;
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSnpPosition() throws IOException
    {
        return this.index >= 0;
    }
    
    /**
     * {@inheritDoc}
     */
    public long getExtentInBasePairs()
    {
        return 1 + this.snps[this.snps.length - 1].getPositionInBasePairs() -
               this.getStartInBasePairs();
    }
    
    /**
     * {@inheritDoc}
     */
    public long getStartInBasePairs()
    {
        return this.snps[0].getPositionInBasePairs();
    }
    
    /**
     * {@inheritDoc}
     */
    public int getChromosomeNumber()
    {
        return this.chromosomeNumber;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return StreamDirection.REVERSE;
    }
}
