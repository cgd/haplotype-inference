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

package org.jax.haplotype.io;

import java.io.IOException;

import org.jax.geneticutil.data.SingleNucleotidePolymorphism;
import org.jax.geneticutil.data.StrainChromosome;

/**
 * A SNP input stream that's a simple wrapper around a pair of
 * {@link StrainChromosome}s. This is the same as
 * {@link ReverseStrainChromosomeSnpInputStream} except we're reading
 * forward.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class ForwardStrainChromosomeSnpInputStream implements SnpInputStream
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
    public ForwardStrainChromosomeSnpInputStream(
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
    public ForwardStrainChromosomeSnpInputStream(
            SingleNucleotidePolymorphism[] referenceSnps,
            SingleNucleotidePolymorphism[] comparisonSnps)
    {
        this.referenceSnps = referenceSnps;
        this.comparisonSnps = comparisonSnps;
        this.index = 0;
    }

    /**
     * {@inheritDoc}
     */
    public boolean getNextSnp() throws IOException
    {
        boolean nextSnp = this.referenceSnps[this.index].getSnpType() ==
                  this.comparisonSnps[this.index].getSnpType();
        this.index++;
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
        return this.index < this.referenceSnps.length;
    }
    
    /**
     * {@inheritDoc}
     */
    public StreamDirection getReadDirection() throws IOException
    {
        return StreamDirection.FORWARD;
    }
}
