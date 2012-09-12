/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.alignment.hmm.tools;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.graph.filter.BloomFilter.RightCodonFacade;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;

/**
 *
 * @author fishjord
 */
public class BloomStats {

    /*
    static String[] tabs;
    static int totalWords;
    static int[] wordsAtDepth;
    static int[] branches = new int[5];
    static int[][] depthBranches;

    private static void processWord(String word, BloomFilter filter, int depth, int maxDepth) {
        if (depth >= maxDepth) {
            return;
        }

        System.out.println(tabs[depth] + word);

        wordsAtDepth[depth]++;
        totalWords++;

        RightCodonFacade walker = filter.new RightCodonFacade(word, 0);
        if (!walker.probeRight()) {
            branches[0]++;
            depthBranches[depth][0]++;
            return;
        }

        int branchCount = 0;
        do {
            String nextWord = word.substring(1) + walker.getPathString();
            processWord(nextWord, filter, depth + 1, maxDepth);
            branchCount++;
        } while (walker.replaceRight());

        branches[branchCount]++;
        depthBranches[depth][branchCount]++;
    }

    public static void main(String[] args) throws Exception {

        if (args.length != 4 && args.length != 3) {
            System.err.println("BloomStats <binary_bloom_filter> <starting word> <max_depth>");
            System.exit(1);
        }

        File inFile = new File(args[0]);
        String startingWord = args[1];
        int wordSize = startingWord.length();
        int maxDepth = Integer.valueOf(args[2]) + 1;

        BloomFilter filter;
        System.err.println("Binary bloom filter file exists, loading from cache");
        long startTime = System.currentTimeMillis();
        ObjectInputStream ois = null;
        ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(inFile)));
        filter = (BloomFilter) ois.readObject();
        ois.close();
        System.err.println("Bloomfilter read from cache in " + (System.currentTimeMillis() - startTime) + " ms");

        System.err.println("Word size= " + wordSize);
        System.err.println("Starting word= " + startingWord);

        wordsAtDepth = new int[maxDepth];
        depthBranches = new int[maxDepth][5];
        tabs = new String[maxDepth];
        for (int index = 0; index < maxDepth; index++) {
            tabs[index] = "";
            for (int i = 0; i < index; i++) {
                tabs[index] += " ";
            }
        }

        processWord(startingWord, filter, 0, maxDepth);

        System.err.println();
        System.err.println("Total words seen: " + totalWords);
        System.err.println();
        System.err.println("Depth\tWords\t0 branches\t1 branches\t2 branches\t3 branches\t4 branches");
        for (int index = 0; index < maxDepth; index++) {
            System.err.println(index + "\t" + wordsAtDepth[index] + "\t" + depthBranches[index][0] + "\t" + depthBranches[index][1] + "\t" + depthBranches[index][2] + "\t" + depthBranches[index][3] + "\t" + depthBranches[index][4]);
        }

        System.err.println();
        System.err.println("Branches\tCounts");
        System.err.println(0 + "\t" + branches[0]);
        System.err.println(1 + "\t" + branches[1]);
        System.err.println(2 + "\t" + branches[2]);
        System.err.println(3 + "\t" + branches[3]);
        System.err.println(4 + "\t" + branches[4]);
    }
     * */
     
}
