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

package org.jax.haplotype.inference;

import java.io.IOException;
import java.util.List;

import org.jax.geneticutil.data.PartitionedInterval;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;

/**
 * Interface that haplotype estimation algorithms should implement
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface HaplotypeEstimator
{
    /**
     * Estimate haplotype blocks from the given position and SDP input streams
     * @param sdpInputStream
     *          the SDPs
     * @param positionInputStream
     *          the physical positions of the SNPs
     * @return
     *          the list of SNP blocks
     * @throws IOException
     *          if we catch an exception reading from the streams
     */
    public abstract List<PartitionedInterval> estimateHaplotypeBlocks(
            SdpInputStream sdpInputStream,
            SnpPositionInputStream positionInputStream) throws IOException;
}
