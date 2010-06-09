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
public interface SnpPositionEvaluator
{
    /**
     * Evaluate the chromosomes at the given snp position
     * @param chromosomesToCompare
     *          the chromosomes we're evaluating
     * @param currSnpIndex
     *          the position that we're evaluating at
     * @return
     *          true if it evaluates as true
     */
    boolean evaluate(StrainChromosome[] chromosomesToCompare, int currSnpIndex);

}
