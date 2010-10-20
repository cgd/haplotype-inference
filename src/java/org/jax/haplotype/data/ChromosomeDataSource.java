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

package org.jax.haplotype.data;

import java.io.Serializable;
import java.util.Set;

import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.SdpInputStream;
import org.jax.haplotype.io.SnpPositionInputStream;
import org.jax.haplotype.io.StreamDirection;

/**
 * The chromosome data source
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public interface ChromosomeDataSource extends Serializable
{
    /**
     * Getter for the genotype data
     * @param strainsToParse
     *          the strains to parse (null means use all available)
     * @return
     *          the data
     */
    public Set<StrainChromosome> getGenotypeData(
            Set<String> strainsToParse);
    
    /**
     * Get the available strains
     * @return
     *          the strain name set
     */
    public Set<String> getAvailableStrains();
    
    /**
     * Get the chromosome number for this data source
     * @return
     *          the chromosome number
     */
    public int getChromosomeNumber();
    
    /**
     * Get the starting position in base pairs
     * @return
     *          the starting position
     */
    public long getDataStartInBasePairs();
    
    /**
     * Getter for the extent
     * @return
     *          the extent in base pairs
     */
    public long getDataExtentInBasePairs();
    
    /**
     * Get SDPs for the strain names. {@link StreamDirection#FORWARD} is
     * assumed
     * @param strainNames
     *          the strain names
     * @return
     *          the SDP stream
     */
    public SdpInputStream getSdpInputStream(String[] strainNames);
    
    /**
     * Get SDPs for the strain names.
     * @param streamDirection
     *          the stream direction to use
     * @param strainNames
     *          the strain names
     * @return
     *          the SDP stream
     */
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String[] strainNames);
    
    /**
     * Get SDPs for the strain names vs the reference strain.
     * {@link StreamDirection#FORWARD} is assumed
     * @param referenceStrainName
     *          the reference strain. the SDPs generated should have a
     *          1 if they are equal to this strain and a 0 otherwise
     * @param comparisonStrainNames
     *          the comparison strain names
     * @return
     *          the SDP stream
     */
    public SdpInputStream getSdpInputStream(
            String referenceStrainName,
            String[] comparisonStrainNames);
    
    /**
     * Get SDPs for the strain names vs the reference strain.
     * {@link StreamDirection#FORWARD} is assumed
     * @param streamDirection
     *          the stream direction to use
     * @param referenceStrainName
     *          the reference strain. the SDPs generated should have a
     *          1 if they are equal to this strain and a 0 otherwise
     * @param comparisonStrainNames
     *          the comparison strain names
     * @return
     *          the SDP stream
     */
    public SdpInputStream getSdpInputStream(
            StreamDirection streamDirection,
            String referenceStrainName,
            String[] comparisonStrainNames);
    
    /**
     * Get the SNP positions stream.
     * {@link StreamDirection#FORWARD} is assumed
     * @return
     *          the SNP positions stream
     */
    public SnpPositionInputStream getSnpPositionInputStream();
    
    /**
     * Get the SNP positions stream.
     * {@link StreamDirection#FORWARD} is assumed
     * @param streamDirection
     *          the stream direction to use
     * @return
     *          the SNP positions stream
     */
    public SnpPositionInputStream getSnpPositionInputStream(
            StreamDirection streamDirection);
}
