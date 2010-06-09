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

package org.jax.haplotype.expressions;

import org.jax.geneticutil.data.StrainChromosome;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class SnpPositionsMatchEvaluator implements SnpPositionEvaluator
{
    /**
     * {@inheritDoc}
     */
    public boolean evaluate(StrainChromosome[] chromosomesToCompare, int snpPosition)
    {
        for(int i = 1; i < chromosomesToCompare.length; i++)
        {
            if(chromosomesToCompare[i].getSingleNucleotidePolymorphisms()[snpPosition].getSnpType() !=
               chromosomesToCompare[i - 1].getSingleNucleotidePolymorphisms()[snpPosition].getSnpType())
            {
                return false;
            }
        }

        return true;
    }
}
