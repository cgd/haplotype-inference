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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;

import org.jax.geneticutil.data.StrainChromosome;
import org.jax.util.TextWrapper;

/**
 * Parses genotype expressions
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class FunctionalExpressionParser
{
    private final Set<StrainChromosome> chromosomes;
    
    private static final char PARAMETERS_START = '(';

    private static final char PARAMETERS_END = ')';
    
    private static final char PARAMETERS_SEPARATOR = ',';
    
    private static final char EQUALITY_SEPARATOR = '=';
    
    /**
     * Different function types that the parser can handle
     */
    public static enum FunctionType
    {
        /**
         * not
         */
        NOT_FUNCTION
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionName()
            {
                return "not";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionDescription()
            {
                return "The \"not\" function should be used " +
                	   "to invert the equality relationship. This function " +
                	   "takes a single chained chromosome comparison as input and " +
                	   "produces a set of SNP groups as output.\n" +
                	   "Eg: not(A=B=C) means all the SNP locations where at " +
                	   "least one of the given chromosomes (A, B or C) has " +
                	   "a different nucleotied value than the others.";
            }
        },
        
        /**
         * union
         */
        UNION_FUNCTION
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionName()
            {
                return "union";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionDescription()
            {
                return "The \"union\" function should " +
                	   "be used to create a union of SNP positions. This " +
                	   "function takes any number of parameters which can be " +
                	   "either chained chromosome comparisons or " +
                	   "SNP groups. This function produces a SNP group as " +
                	   "output\n" +
                	   "Eg: union(A=B, not(C=B)) will select all of the SNP " +
                	   "positions where either A is the same as B, or C is different " +
                	   "from B.";
            }
        },
        
        /**
         * intersection
         */
        INTERSECTION_FUNCTION
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionName()
            {
                return "intersection";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getFunctionDescription()
            {
                return  "The \"intersection\" function should " +
                        "be used to create an intersection of SNP positions. This " +
                        "function takes any number of parameters which can be " +
                        "either chained chromosome comparisons or " +
                        "SNP groups. This function produces a SNP group as " +
                        "output\n" +
                        "Eg: intersection(A=B, not(C=B)) will select all of the SNP " +
                        "positions where both A is the same as B, and C is different " +
                        "from B.";
            }
        };
        
        /**
         * The name that should be used for the given function
         * @return
         *          the function name
         */
        public abstract String getFunctionName();
        
        /**
         * Get the text description for this function.
         * @return
         *          the description
         */
        public abstract String getFunctionDescription();
        
        /**
         * {@inheritDoc}
         */
        @Override
        public String toString()
        {
            return this.getFunctionName();
        }
        
        /**
         * Determine if the given function name is a case-insensitive match
         * for this function
         * @param functionName
         *          the function name to test
         * @return
         *          true if the function name is a match
         */
        public boolean matches(String functionName)
        {
            return this.getFunctionName().equalsIgnoreCase(functionName);
        }
    }

    /**
     * @param chromosomes
     */
    public FunctionalExpressionParser(Set<StrainChromosome> chromosomes)
    {
        this.chromosomes = chromosomes;
    }

    /**
     * Parses the given expression
     * @param expressionToParse
     *          the expression string to parse
     * @return
     *          the result
     * @throws MalformedExpressionException
     *          if the expression contains errors
     */
    public SnpIntervalExpression parseExpression(String expressionToParse)
    throws MalformedExpressionException
    {
        int paramStartIndex = expressionToParse.indexOf(PARAMETERS_START);
        
        if(paramStartIndex == -1)
        {
            StrainChromosome[] equalityChromosomes =
                this.parseExpressionAsEqualityComparison(expressionToParse);
            return new SnpPositionExpression(
                    equalityChromosomes,
                    new SnpPositionsMatchEvaluator());
        }
        else
        {
            // break out the function name & parameters
            String functionName =
                expressionToParse.substring(0, paramStartIndex).trim();
            String parametersString =
                expressionToParse.substring(paramStartIndex).trim();
            String[] parameters = this.breakUpParameters(parametersString);
            for(int i = 0; i < parameters.length; i++)
            {
                parameters[i] = parameters[i].trim();
            }
            
            return this.createFunctionExpression(
                    functionName,
                    parameters);
        }
    }

    /**
     * Create a function expression using the function name and parameters
     * @param functionName
     *          the function name to build an expression from
     * @param parameters
     *          the parameters to use
     * @return
     *          the result
     * @throws MalformedExpressionException 
     *          if the expression is bad
     */
    private SnpIntervalExpression createFunctionExpression(
            String functionName,
            String[] parameters)
    throws
            MalformedExpressionException
    {
        if(functionName.length() == 0)
        {
            throw new MalformedExpressionException(
                    "opening braces \"(\" must be proceeded by " +
                    "a function name!");
        }
        else
        {
            if(FunctionType.NOT_FUNCTION.matches(functionName))
            {
                return this.createNotFunctionFromParameters(parameters);
            }
            else
            {
                // convert the arguments into expressions
                SnpIntervalExpression[] parameterExpressions = new SnpIntervalExpression[parameters.length];
                for(int i = 0; i < parameterExpressions.length; i++)
                {
                    parameterExpressions[i] = this.parseExpression(parameters[i]);
                }

                if(FunctionType.INTERSECTION_FUNCTION.matches(functionName))
                {
                    return this.createIntersectionFunctionFromParameters(parameterExpressions);
                }
                else if(FunctionType.UNION_FUNCTION.matches(functionName))
                {
                    return this.createUnionFunctionFromParameters(parameterExpressions);
                }
                else
                {
                    StringBuffer errorMessage = new StringBuffer(
                            "found bad function name: \"" + functionName +
                            "\". Permissable values are:");
                    for(FunctionType currFunction: FunctionType.values())
                    {
                        errorMessage.append(" ");
                        errorMessage.append(currFunction.getFunctionName());
                    }
                    
                    throw new MalformedExpressionException(errorMessage.toString());
                }
            }
        }
    }

    private UnionIntervalExpression createUnionFunctionFromParameters(
            SnpIntervalExpression[] parameterExpressions) throws MalformedExpressionException
    {
        if(parameterExpressions.length < 2)
        {
            throw new MalformedExpressionException(
                    "union must have at least two parameters");
        }
        else
        {
            UnionIntervalExpression unionExpression = new UnionIntervalExpression(
                    parameterExpressions[0],
                    parameterExpressions[1]);
            
            for(int i = 2; i < parameterExpressions.length; i++)
            {
                unionExpression = new UnionIntervalExpression(
                        unionExpression,
                        parameterExpressions[i]);
            }
            
            return unionExpression;
        }
    }

    private IntersectionIntervalExpression createIntersectionFunctionFromParameters(
            SnpIntervalExpression[] parameterExpressions) throws MalformedExpressionException
    {
        if(parameterExpressions.length < 2)
        {
            throw new MalformedExpressionException(
                    "intersection must have at least two parameters");
        }
        else
        {
            IntersectionIntervalExpression intersectionExpression = new IntersectionIntervalExpression(
                    parameterExpressions[0],
                    parameterExpressions[1]);
            
            for(int i = 2; i < parameterExpressions.length; i++)
            {
                intersectionExpression = new IntersectionIntervalExpression(
                        intersectionExpression,
                        parameterExpressions[i]);
            }
            
            return intersectionExpression;
        }
    }

    private SnpIntervalExpression createNotFunctionFromParameters(
            String[] parameters) throws MalformedExpressionException
    {
        if(parameters.length != 1)
        {
            throw new MalformedExpressionException(
                    "the \"not\" function is expected to take a single " +
                    "chromosome equality expression as a parameter");
        }
        else
        {
            StrainChromosome[] equalityChromosomes = this.parseExpressionAsEqualityComparison(
                    parameters[0]);
            NotSnpPositionEvaluator notEvaluator = new NotSnpPositionEvaluator(
                    new SnpPositionsMatchEvaluator());
            
            return new SnpPositionExpression(
                    equalityChromosomes,
                    notEvaluator);
        }
    }

    private String[] breakUpParameters(String parametersString) throws MalformedExpressionException
    {
        if(!parametersString.startsWith("" + PARAMETERS_START))
        {
            throw new MalformedExpressionException(
                    "parameters: " + parametersString +
                    " should start with a " + PARAMETERS_START);
        }
        else if(!parametersString.endsWith("" + PARAMETERS_END))
        {
            throw new MalformedExpressionException(
                    "parameters: " + parametersString +
                    " should end with a " + PARAMETERS_END);
        }
        else
        {
            // trim off the braces
            parametersString = parametersString.substring(1, parametersString.length() - 1);
            
            List<String> parameterStringList = new ArrayList<String>();
            int parameterStartCursor = 0;
            int parameterEndCursor = this.getNextParameterEnd(
                    parametersString,
                    parameterStartCursor);
            while(parameterStartCursor < parameterEndCursor - 1)
            {
                parameterStringList.add(
                        parametersString.substring(
                                parameterStartCursor,
                                parameterEndCursor));
                
                parameterStartCursor = parameterEndCursor + 1;
                parameterEndCursor = this.getNextParameterEnd(
                        parametersString,
                        parameterStartCursor + 1);
            }
            
            return parameterStringList.toArray(
                    new String[parameterStringList.size()]);
        }
    }

    private int getNextParameterEnd(String parametersString,
            int parameterStartCursor) throws MalformedExpressionException
    {
        int startEndParamBalance = 0;
        int parameterEndCursor = parameterStartCursor;
        while(parameterEndCursor < parametersString.length())
        {
            char currChar = parametersString.charAt(parameterEndCursor);
            if(currChar == PARAMETERS_START)
            {
                startEndParamBalance++;
            }
            else if(currChar == PARAMETERS_END)
            {
                startEndParamBalance--;
                if(startEndParamBalance < 0)
                {
                    return parameterEndCursor;
                }
            }
            else
            {
                if(startEndParamBalance == 0 && currChar == PARAMETERS_SEPARATOR)
                {
                    return parameterEndCursor;
                }
            }
            
            parameterEndCursor++;
        }
        
        if(startEndParamBalance != 0)
        {
            throw new MalformedExpressionException(
                    "unmatched braces found in " + parametersString);
        }
        else
        {
            return parametersString.length();
        }
    }

    private StrainChromosome[] parseExpressionAsEqualityComparison(String currLine) throws MalformedExpressionException
    {
        StringTokenizer equalityTokenizer = new StringTokenizer(
                currLine,
                Character.toString(EQUALITY_SEPARATOR));
        
        if(equalityTokenizer.countTokens() < 1)
        {
            throw new MalformedExpressionException(
                    currLine + " is not a valid chromosome equality comparison");
        }
        else
        {
            StrainChromosome[] selectedChromosomes =
                new StrainChromosome[equalityTokenizer.countTokens()];
            for(int i = 0; i < selectedChromosomes.length; i++)
            {
                String currChromosomeName = equalityTokenizer.nextToken();
                StrainChromosome currChromosome = this.getChromosomeByName(currChromosomeName);
                
                if(currChromosome == null)
                {
                    throw new MalformedExpressionException(
                            "could not find strain named: " + currChromosomeName);
                }
                else
                {
                    selectedChromosomes[i] = currChromosome;
                }
            }
            
            return selectedChromosomes;
        }
    }

    private StrainChromosome getChromosomeByName(String currChromosomeName)
    {
        for(StrainChromosome currChromosome: this.chromosomes)
        {
            if(currChromosome.getStrainName().equals(currChromosomeName))
            {
                return currChromosome;
            }
        }
        
        return null;
    }

    /**
     * Print help info for expressions using the wrap margin
     * @param wrapMargin
     *          the margin to use
     */
    public static void printExpressionHelp(int wrapMargin)
    {
        final String prefix = "   ";
        final String msg1 =
            "SNP expressions are used to operate on " +
            "SNP data. The expressions are built with chained chromosome " +
            "comparisons which produce SNP groups as output and functions " +
            "which take the SNP groups as input and transform them in some " +
            "way.";
        FunctionalExpressionParser.printWrapped(msg1, wrapMargin, "");
        final String msg2 =
            "Eg: \"C57BL/6J=C3H/HeJ\" can be used to query for where " +
            "Black-6 is the same as C3H. From this you can build more " +
            "complicated expressions like:";
        FunctionalExpressionParser.printWrapped(msg2, wrapMargin, "");
        System.out.println(
                "intersection(C57BL/6J=C3H/HeJ, union(not(SEG/Pas=SJL/J=SM/J), NZB/BlNJ=NZW/LacJ))");
        System.out.println();
        for(FunctionType currFunction: FunctionType.values())
        {
            System.out.println("Function: " + currFunction.getFunctionName());
            FunctionalExpressionParser.printWrapped(
                    currFunction.getFunctionDescription(),
                    wrapMargin,
                    prefix);
        }
    }
    
    private static void printWrapped(String message, int wrapMargin, String prefix)
    {
        for(String currLine: TextWrapper.wrapText(message, wrapMargin - prefix.length()))
        {
            System.out.println(prefix + currLine);
        }
    }

    /**
     * Getter for the chromosomes that this expression parser operates on
     * @return the chromosomes
     */
    public Set<StrainChromosome> getChromosomes()
    {
        return this.chromosomes;
    }
}
