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
 * A SNP input stream that's a simple wrapper around a pair of
 * {@link StrainChromosome}s. This is the same as
 * {@link ForwardStrainChromosomeSnpInputStream} except we're reading
 * in reverse.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ReverseStrainChromosomeSnpInputStream implements SnpInputStream
{
    private final SingleNucleotidePolymorphism[] referenceSnps;
    
    private final SingleNucleotidePolymorphism[] comparisonSnps;
    
    private int index;
    
    /**
     * Constructor
     * @param referenceStrainChromosome
     *          the reference strain
     * @param comparisonStrainChromosome
     *          the comparison strain
     */
    public ReverseStrainChromosomeSnpInputStream(
            StrainChromosome referenceStrainChromosome,
            StrainChromosome comparisonStrainChromosome)
    {
        this(referenceStrainChromosome.getSingleNucleotidePolymorphisms(),
             comparisonStrainChromosome.getSingleNucleotidePolymorphisms());
    }

    /**
     * Constructor
     * @param referenceSnps
     *          the reference snp values
     * @param comparisonSnps
     *          the snp values that we're comparing to the reference
     */
    public ReverseStrainChromosomeSnpInputStream(
            SingleNucleotidePolymorphism[] referenceSnps,
            SingleNucleotidePolymorphism[] comparisonSnps)
    {
        this.referenceSnps = referenceSnps;
        this.comparisonSnps = comparisonSnps;
        this.index = referenceSnps.length - 1;
    }

    /**
     * {@inheritDoc}
     */
    public boolean getNextSnp() throws IOException
    {
        boolean nextSnp =
                  this.referenceSnps[this.index].getSnpType() ==
                  this.comparisonSnps[this.index].getSnpType();
        this.index--;
        return nextSnp;
    }

    /**
     * {@inheritDoc}
     */
    public long getSnpCount() throws IOException
    {
        return this.referenceSnps.length;
    }

    /**
     * {@inheritDoc}
     */
    public boolean hasNextSnp() throws IOException
    {
        return this.index >= 0;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return StreamDirection.REVERSE;
    }
}
