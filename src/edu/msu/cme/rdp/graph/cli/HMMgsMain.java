/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.cli;

import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class HMMgsMain {

    private static void printUsageAndExit() {
        System.err.println("USAGE: HMMgs <command> <options>");
        System.err.println("\tbuild       - Build a bloom filter");
        System.err.println("\tstats       - Display bloom filter stats");
        System.err.println("\tsearch      - Search a bloom filter with an hmm");
        System.err.println("\tmerge       - Merge HMMgs left and right fragments");
        System.exit(1);
    }

    public static void main(String[] args) throws Exception {
        if(args.length == 0) {
            printUsageAndExit();
        }

        String cmd = args[0];
        args = Arrays.copyOfRange(args, 1, args.length);

        if(cmd.equals("build")) {
            BloomFilterBuilder.main(args);
        } else if(cmd.equals("append")) {
            BloomFilterAppender.main(args);
        } else if (cmd.equals("stats")) {
            BloomFilterStats.main(args);
        } else if(cmd.equals("search")) {
            TimeLimitedSearch.main(args);
        } else if(cmd.equals("merge")) {
            ContigMerger.main(args);
        } else {
            printUsageAndExit();
        }
    }
}
