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
