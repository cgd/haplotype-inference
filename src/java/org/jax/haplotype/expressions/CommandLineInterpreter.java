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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.Set;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.StringTokenizer;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import javax.swing.JFileChooser;

import org.jax.geneticutil.data.BasePairInterval;
import org.jax.geneticutil.data.StrainChromosome;
import org.jax.haplotype.io.GenotypeParser;
import org.jax.util.TextWrapper;
import org.jax.util.TypeSafeSystemProperties;

/**
 * A command line interpreter for SNP intervals.
 * @author <A HREF="mailto:keith.sheppard@jax.org">Keith Sheppard</A>
 */
public class CommandLineInterpreter
{
    /**
     * our logger
     */
    private static final Logger LOG = Logger.getLogger(
            CommandLineInterpreter.class.getName());
    
    private static final int LINE_WRAP_COLUMN = 120;
    
    /**
     * The command enumeration
     */
    private static enum Command
    {
        /**
         * help command
         */
        HELP
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandHelp()
            {
                return "Print this help.";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandName()
            {
                return "help";
            }
        },
        
        /**
         * expression help command
         */
        EXPRESSION_HELP
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandHelp()
            {
                return "Print help information for expression the expression " +
                	   "syntax.";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandName()
            {
                return "expression-help";
            }
        },
        
        /**
         * show the strains
         */
        STRAINS
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandHelp()
            {
                return "List the names of all loaded strains.";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandName()
            {
                return "strains";
            }
        },
        
        /**
         * command for loading genotype data
         */
        LOAD_DATA
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandHelp()
            {
                return "Load strain data. Currently the only supported data " +
                	   "format is a specific comma separated value format. " +
                	   "To see more information, issue the \"" +
                	   this.getCommandName() + "\" command.";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandName()
            {
                return "load";
            }
        },
        
        /**
         * quit the interpreter
         */
        QUIT
        {
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandHelp()
            {
                return "Quit the SNP expression interpreter.";
            }
            
            /**
             * {@inheritDoc}
             */
            @Override
            public String getCommandName()
            {
                return "quit";
            }
        };
        
        /**
         * Get the command name
         * @return
         *          the command name
         */
        public abstract String getCommandName();
        
        /**
         * Get the text for this command's help
         * @return
         *          the text
         */
        public abstract String getCommandHelp();
        
        /**
         * Print help for this command
         */
        public void printHelp()
        {
            System.out.println(this.getCommandName() + ":");
            final String prefix = "   ";
            String[] wrappedText =
                TextWrapper.wrapText(
                        this.getCommandHelp(),
                        LINE_WRAP_COLUMN - prefix.length());
            for(String currLine: wrappedText)
            {
                System.out.println(prefix + currLine);
            }
        }
        
        /**
         * Print help info for all of the commands.
         */
        public static void printAllHelp()
        {
            for(Command currCommand: Command.values())
            {
                currCommand.printHelp();
            }
        }
        
        /**
         * Get the command with the given name
         * @param commandName
         *          the command name
         * @return
         *          the matching command
         */
        public static Command getCommandWithName(String commandName)
        {
            for(Command currCommand: Command.values())
            {
                if(currCommand.getCommandName().equals(commandName))
                {
                    return currCommand;
                }
            }
            
            return null;
        }
    }
    
    /**
     * @param dataFileString
     * @param expressionString
     * @param outputFileString
     */
    private void runInterpreter(String dataFileString, String expressionString, String outputFileString)
    {
        try
        {
            GenotypeParser dataFileParser = new GenotypeParser();
            Set<StrainChromosome> genotypeData =
                dataFileParser.parseGenotypeFromStream(new FileInputStream(dataFileString));
            
            this.parseExpression(
                    expressionString,
                    new FunctionalExpressionParser(genotypeData),
                    outputFileString);
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "failed to parse expression",
                    ex);
        }
    }
    
    /**
     * Run the command line interpreter.
     */
    public void runInterpreter()
    {
        this.runInterpreter(null);
    }
    
    private String readDataFileName() throws IOException
    {
        System.out.println(
                "Please enter the name of a raw SNP data file. The following " +
                "gives an example of what an example header and some data " +
        "rows for this file should look like:");
        System.out.println();
        System.out.println(
        "snp.ID,ChrID,build.36.bp.Position,Source,129S1/SvImJ,129S4/SvJae,129X1/SvJ");
        System.out.println(
        "NES15061566,chr10,3011171,Perlegen36_b03,t,c,c");
        System.out.println(
        "rs29367782|NES15061567,chr10,3011241,Celera2-Perlegen36_b03,A,c,C");
        System.out.println();
        return this.readFileString();
    }
    
    /**
     * Run the command line interpreter using the given data file.
     * @param dataFileString
     */
    public void runInterpreter(String dataFileString)
    {
        GenotypeParser dataFileParser = new GenotypeParser();
        
        Set<StrainChromosome> genotypeData = null;
        FunctionalExpressionParser parser = null;
        
        try
        {
            if(dataFileString != null)
            {
                genotypeData =
                    dataFileParser.parseGenotypeFromStream(new FileInputStream(dataFileString));
                parser = new FunctionalExpressionParser(genotypeData);
            }
            
            // Read next line until we hit an end-of-file or we're killed
            for(String currLine = this.readExpressionString();
                currLine != null;
                currLine = this.readExpressionString())
            {
                // see if it's a command before we treat it as an
                // expression
                StringTokenizer currLineTok = new StringTokenizer(currLine);
                if(currLineTok.countTokens() > 0)
                {
                    String potentialCommand = currLineTok.nextToken();
                    Command matchingCommand =
                        Command.getCommandWithName(potentialCommand);
                    if(matchingCommand == null)
                    {
                        // default behavior is to treat it as an
                        // expression
                        if(parser == null)
                        {
                            this.sayNoGenotypeDataIsLoaded();
                        }
                        else
                        {
                            this.parseExpression(currLine, parser);
                        }
                    }
                    else
                    {
                        switch(matchingCommand)
                        {
                            case EXPRESSION_HELP:
                            {
                                FunctionalExpressionParser.printExpressionHelp(LINE_WRAP_COLUMN);
                            }
                            break;
                            
                            case HELP:
                            {
                                Command.printAllHelp();
                            }
                            break;
                            
                            case LOAD_DATA:
                            {
                                dataFileString = this.readDataFileName();
                                
                                genotypeData = null;
                                parser = null;
                                
                                if(dataFileString != null)
                                {
                                    genotypeData =
                                        dataFileParser.parseGenotypeFromStream(new FileInputStream(dataFileString));
                                    parser = new FunctionalExpressionParser(genotypeData);
                                }
                            }
                            break;
                            
                            case QUIT:
                            {
                                // just return
                                return;
                            }
                            
                            case STRAINS:
                            {
                                if(genotypeData == null)
                                {
                                    this.sayNoGenotypeDataIsLoaded();
                                }
                                else
                                {
                                    SortedMap<String, SortedSet<Integer>> strainToChromoMap =
                                        new TreeMap<String, SortedSet<Integer>>();
                                    for(StrainChromosome chromosome: genotypeData)
                                    {
                                        SortedSet<Integer> chromoNums =
                                            strainToChromoMap.get(chromosome.getStrainName());
                                        if(chromoNums == null)
                                        {
                                            chromoNums = new TreeSet<Integer>();
                                            strainToChromoMap.put(
                                                    chromosome.getStrainName(),
                                                    chromoNums);
                                        }
                                        
                                        chromoNums.add(chromosome.getChromosomeNumber());
                                    }
                                    
                                    for(Entry<String, SortedSet<Integer>> chromoEntry: strainToChromoMap.entrySet())
                                    {
                                        System.out.println(chromoEntry.getKey());
                                        System.out.print("   Chromosome Numbers: ");
                                        for(Integer currChromoNum: chromoEntry.getValue())
                                        {
                                            System.out.print(currChromoNum.intValue());
                                        }
                                        System.out.println();
                                    }
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "error occured while processing commands",
                    ex);
        }
    }
    
    private void sayNoGenotypeDataIsLoaded()
    {
        System.out.println(
                "There is no genotype data loaded. Use the \"" +
                Command.LOAD_DATA.getCommandName() +
                "\" command.");
    }
    
    private void parseExpression(
            String expressionString,
            FunctionalExpressionParser parser)
    {
        this.parseExpression(expressionString, parser, null);
    }
    
    private void parseExpression(
            String expressionString,
            FunctionalExpressionParser parser,
            String saveFileString)
    {
        try
        {
            // Convert the current line into an expression
            SnpIntervalExpression currExpression = parser.parseExpression(expressionString);
            BasePairInterval[] snpIntervals =
                currExpression.evaluateExpression();
            
            // save the results to file
            if(saveFileString == null)
            {
                System.out.println(
                        "The expression has been parsed and has " +
                        "detected " + snpIntervals.length +
                        " matching SNP intervals. Please enter a file to save " +
                        "the intervals as a comma separated file: ");
                saveFileString = this.readFileString();
            }
            
            if(saveFileString == null)
            {
                System.out.println("Discarding results of previous expression.");
            }
            else
            {
                try
                {
                    File saveFile = new File(saveFileString);
                    saveFile.createNewFile();
                    PrintStream saveFileStream = new PrintStream(saveFile);
                    saveFileStream.println("# Expression Parsed: " + expressionString);
                    saveFileStream.println("# Loaded Chromosomes:");
                    Set<StrainChromosome> chromosomes = parser.getChromosomes();
                    for(StrainChromosome currChromosomes: chromosomes)
                    {
                        saveFileStream.println(
                                "#   Strain: " + currChromosomes.getStrainName() +
                                " Chromosome #: " + currChromosomes.getChromosomeNumber());
                    }
                    saveFileStream.println("#");
                    saveFileStream.println("# snp starting index (0 based), extent in snps");
                    for(BasePairInterval snpInterval: snpIntervals)
                    {
                        saveFileStream.println(
                                snpInterval.getStartInBasePairs() +
                                "," +
                                snpInterval.getExtentInBasePairs());
                    }
                    saveFileStream.flush();
                    saveFileStream.close();
                }
                catch(Exception ex)
                {
                    LOG.log(Level.SEVERE,
                            "Failed to save expression results to: " + saveFileString,
                            ex);
                }
            }
        }
        catch(Exception ex)
        {
            LOG.log(Level.SEVERE,
                    "Failed to parse: " + expressionString,
                    ex);
        }
    }
    
    /**
     * Get the next expression string from user input.
     * @return  the next expression string
     * @throws IOException  if stdin dies on us
     */
    private String readExpressionString() throws IOException
    {
        BufferedReader bufferedIn = new BufferedReader(
                new InputStreamReader(System.in));
        String expressionString;
        
        do
        {
            // Print the prompt and read in the expression
            // (if the expression is empty, keep trying)
            System.out.print("> ");
            expressionString = bufferedIn.readLine();
            if(expressionString != null)
            {
                expressionString = expressionString.trim();
            }
        } while(expressionString != null && expressionString.length() == 0);
        
        if(expressionString == null || expressionString.equalsIgnoreCase("quit"))
        {
            return null;
        }
        else
        {
            return expressionString;
        }
    }
    
    /**
     * Get the next expression string from user input.
     * @return  the next expression string
     * @throws IOException  if stdin dies on us
     */
    private String readFileString() throws IOException
    {
        BufferedReader bufferedIn = new BufferedReader(
                new InputStreamReader(System.in));
        String fileString;
        
        do
        {
            // Print the prompt and read in the expression
            // (if the expression is empty, keep trying)
            System.out.print("file name, \"browse\" or \"cancel\"> ");
            fileString = bufferedIn.readLine();
            if(fileString != null)
            {
                fileString = fileString.trim();
                if(fileString.equalsIgnoreCase("browse"))
                {
                    JFileChooser fileChooser = new JFileChooser();
                    fileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
                    fileChooser.setMultiSelectionEnabled(false);
                    fileChooser.setCurrentDirectory(
                            new File(TypeSafeSystemProperties.getWorkingDirectory()));
                    int userSelection = fileChooser.showOpenDialog(null);
                    if(userSelection == JFileChooser.APPROVE_OPTION)
                    {
                        File selectedFile = fileChooser.getSelectedFile();
                        fileString = selectedFile.getAbsolutePath();
                        System.out.println(fileString);
                    }
                    else
                    {
                        // loop again
                        fileString = "";
                        System.out.println("browse canceled");
                    }
                }
            }
        } while(fileString != null && fileString.length() == 0);
        
        if(fileString == null || fileString.equalsIgnoreCase("cancel"))
        {
            return null;
        }
        else
        {
            return fileString;
        }
    }
    
    /**
     * The main function to start our interpreter.
     * @param args  the application arguments (not used)
     */
    public static void main(String[] args)
    {
        CommandLineInterpreter cli = new CommandLineInterpreter();
        if(args.length == 1)
        {
            cli.runInterpreter(args[0]);
        }
        else if(args.length == 3)
        {
            cli.runInterpreter(args[0], args[1], args[2]);
        }
        else
        {
            cli.runInterpreter();
        }
    }

}
