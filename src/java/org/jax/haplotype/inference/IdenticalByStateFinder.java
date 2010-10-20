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

package org.jax.haplotype.inference;

import java.util.Set;

import org.jax.geneticutil.data.SnpIntervalList;
import org.jax.geneticutil.data.SnpIntervalListGroup;
import org.jax.geneticutil.data.StrainChromosome;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface IdenticalByStateFinder
{
    /**
     * Find all of the regions that are IBS
     * @param chromosome1
     *          the 1st chromosome
     * @param chromosome2
     *          the 2nd chromosome
     * @param minimumExtentInSnps
     *          the minimum extent in snps
     * @param minimumExtentInBasePairs
     *          the minimum extent in base pairs
     * @return
     *          the result
     */
    public SnpIntervalList findIdenticalByStateRegions(
            StrainChromosome chromosome1,
            StrainChromosome chromosome2,
            int minimumExtentInSnps,
            long minimumExtentInBasePairs);
    
    /**
     * Find all of the regions that are IBS
     * @param referenceStrainChromosome
     *          the reference strain chromosome
     * @param comparisonStrainChromosomes
     *          the comparison strain chromosomes
     * @param minimumExtentInSnps
     *          the minimum extent in snps
     * @param minimumExtentInBasePairs
     *          the minimum extent in base pairs
     * @return
     *          the result
     */
    public SnpIntervalListGroup findIdenticalByStateRegions(
            final StrainChromosome referenceStrainChromosome,
            final Set<StrainChromosome> comparisonStrainChromosomes,
            final int minimumExtentInSnps,
            final long minimumExtentInBasePairs);
}
