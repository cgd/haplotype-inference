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

import org.jax.haplotype.jaxbgenerated.ChromosomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.CommaSeparatedChromosomeDataSourceType;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class JaxbChromosomeDataSourceFactory
{
    /**
     * Create a chromosome data source from the given jaxb definition
     * @param jaxbChromosomeDataSource
     *          the jaxb definition
     * @return
     *          the data source
     */
    public static final ChromosomeDataSource getChromosomeDataSource(
            ChromosomeDataSourceType jaxbChromosomeDataSource)
    {
        if(jaxbChromosomeDataSource instanceof CommaSeparatedChromosomeDataSourceType)
        {
            CommaSeparatedChromosomeDataSourceType jaxbCsvChromoDataSourceEntity =
                (CommaSeparatedChromosomeDataSourceType)jaxbChromosomeDataSource;
            return new CommaSeparatedChromosomeDataSource(
                    new File(jaxbCsvChromoDataSourceEntity.getFileLocation()),
                    jaxbCsvChromoDataSourceEntity.getChromosomeNumber());
        }
        else
        {
            throw new RuntimeException(
                    "Unknown chromosome data source type: " +
                    jaxbChromosomeDataSource.getClass().getName());
        }
    }
}
