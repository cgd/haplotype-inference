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

import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.xml.bind.JAXBContext;

import org.jax.haplotype.data.GenomeDataSource;
import org.jax.haplotype.data.JaxbGenomeDataSourceFactory;
import org.jax.haplotype.jaxbgenerated.GenomeDataSourceType;
import org.jax.haplotype.jaxbgenerated.GenomeDataSources;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CachingGenomeDataManager
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            CachingGenomeDataManager.class.getName());
    
    private static final CachingGenomeDataManager instance =
        new CachingGenomeDataManager();
    
    private static final String GENOME_METADATA_RESOURCE_LOCATION =
        "/genome-data-sources.xml";
    
    private final Map<String, GenomeDataSource> genomeDataMap =
        new HashMap<String, GenomeDataSource>();
    
    /**
     * 
     */
    private CachingGenomeDataManager()
    {
        try
        {
            JAXBContext jaxbContext = JAXBContext.newInstance(
                    GenomeDataSources.class.getPackage().getName());
            InputStream genomeMetadataStream = CachingGenomeDataManager.class.getResourceAsStream(
                    GENOME_METADATA_RESOURCE_LOCATION);
            GenomeDataSources jaxbGenomeDataSources =
                (GenomeDataSources)jaxbContext.createUnmarshaller().unmarshal(
                        genomeMetadataStream);
            
            for(GenomeDataSourceType jaxbGenomeDataSource:
                jaxbGenomeDataSources.getGenomeDataSource())
            {
                GenomeDataSource genomeDataSource =
                    JaxbGenomeDataSourceFactory.getGenomeDataSource(jaxbGenomeDataSource);
                this.genomeDataMap.put(
                        genomeDataSource.getName(),
                        genomeDataSource);
            }
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to read jaxb data",
                    ex);
        }
    }
    
    /**
     * Getter for the instance
     * @return
     *          the instance
     */
    public static CachingGenomeDataManager getInstance()
    {
        return CachingGenomeDataManager.instance;
    }
    
    /**
     * Getter for the genome data map where keys are genome names
     * @return
     *          the genome data map
     */
    public Map<String, GenomeDataSource> getGenomeDataMap()
    {
        return this.genomeDataMap;
    }
}
