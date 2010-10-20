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

package org.jax.haplotype.expressions;

/**
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class MalformedExpressionException extends Exception
{
    /**
     * every {@link java.io.Serializable} should have one of these
     */
    private static final long serialVersionUID = 6446308254056374233L;

    /**
     * 
     */
    public MalformedExpressionException()
    {
        super();
    }

    /**
     * @param message
     * @param cause
     */
    public MalformedExpressionException(String message, Throwable cause)
    {
        super(message, cause);
    }

    /**
     * @param message
     */
    public MalformedExpressionException(String message)
    {
        super(message);
    }

    /**
     * @param cause
     */
    public MalformedExpressionException(Throwable cause)
    {
        super(cause);
    }
}
