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

package org.jax.haplotype.phylogeny.data;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class NoValidPhylogenyException extends Exception
{
    /**
     * every {@link java.io.Serializable} is supposed to have one of these
     */
    private static final long serialVersionUID = -8310191072630954004L;

    /**
     * Constructor
     */
    public NoValidPhylogenyException()
    {
        super();
    }

    /**
     * Constructor
     * @param message
     *          the error message
     */
    public NoValidPhylogenyException(String message)
    {
        super(message);
    }

    /**
     * Constructor
     * @param cause
     *          the root cause of this exception
     */
    public NoValidPhylogenyException(Throwable cause)
    {
        super(cause);
    }

    /**
     * Constructor
     * @param message
     *          the error message
     * @param cause
     *          the root cause of this exception
     */
    public NoValidPhylogenyException(String message, Throwable cause)
    {
        super(message, cause);
    }
}
