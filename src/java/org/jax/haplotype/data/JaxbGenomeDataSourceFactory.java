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

package org.jax.haplotype.data;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

import org.jax.haplotype.io.SnpStreamUtil;
import org.jax.haplotype.jaxbgenerated.BinaryGenomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.ChromosomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.GenomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.SimpleGenomeDataSourceType;

/**
 * A factory for JAXB defined genome data sources
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbGenomeDataSourceFactory
{
    /**
     * factory method for the genome data source
     * @param genomeDataSourceElement
     *          the JAXB definition for the genome data
     * @return
     *          the genome data source
     */
    public static final GenomeDataSource getGenomeDataSource(
            GenomeDataSourceType genomeDataSourceElement)
    {
        if(genomeDataSourceElement instanceof SimpleGenomeDataSourceType)
        {
            SimpleGenomeDataSourceType simpleGenomeDataSourceElement =
                (SimpleGenomeDataSourceType)genomeDataSourceElement;
            Map<Integer, ChromosomeDataSource> chromosomeDataSourceMap =
                new HashMap<Integer, ChromosomeDataSource>();
            
            for(ChromosomeDataSourceType jaxbChromosomeDataSource:
                simpleGenomeDataSourceElement.getChromosomeDataSource())
            {
                ChromosomeDataSource chromosomeDataSource =
                    JaxbChromosomeDataSourceFactory.getChromosomeDataSource(
                            jaxbChromosomeDataSource);
                chromosomeDataSourceMap.put(
                        chromosomeDataSource.getChromosomeNumber(),
                        chromosomeDataSource);
            }
            
            return new GenomeDataSource(
                    simpleGenomeDataSourceElement.getName(),
                    simpleGenomeDataSourceElement.getNcbiBuildVersion(),
                    chromosomeDataSourceMap);
        }
        else if(genomeDataSourceElement instanceof BinaryGenomeDataSourceType)
        {
            BinaryGenomeDataSourceType binaryGenomeDataSourceElement =
                (BinaryGenomeDataSourceType)genomeDataSourceElement;
            
            File dataDirectory = new File(
                    binaryGenomeDataSourceElement.getDataDirectory());
            Map<Integer, StreamingBinaryChromosomeDataSource> chromosomeDataSourceMap =
                SnpStreamUtil.getBinaryChromosomeDataSources(
                        dataDirectory);
            
            return new GenomeDataSource(
                    binaryGenomeDataSourceElement.getName(),
                    binaryGenomeDataSourceElement.getNcbiBuildVersion(),
                    chromosomeDataSourceMap);
        }
        else
        {
            throw new RuntimeException(
                    "Unknown genome data source type: " +
                    genomeDataSourceElement.getClass().getName());
        }
    }
}
