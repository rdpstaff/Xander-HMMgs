/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.cli;

import edu.msu.cme.rdp.graph.search.SearchResult;
import java.io.PrintStream;
import java.text.DecimalFormat;

/**
 *
 * @author fishjord
 */
public class HMMBloomSearch {

    private static final DecimalFormat format = new DecimalFormat("#.##");
    private static final DecimalFormat threeDig = new DecimalFormat("#.###");

    public static void printHeader(PrintStream out) {
        printHeader(out, false);
    }

    public static void printHeader(PrintStream out, boolean protSearch) {
        out.println("#contig_id\tstarting kmer\tstarting state\tnucl length" + (protSearch ? "\tprot length" : "") + "\tsearch direction\tnats\tbits\ttime (s)");
    }

    public static void printResult(String seqid, SearchResult result, PrintStream out) {
        printResult(seqid, true, result, out);
    }

    public static void printResult(String seqid, boolean protSearch, SearchResult result, PrintStream out) {
        System.out.println(
                seqid + "\t"
                + result.getKmer() + "\t"
		+ result.getModelPosition() + "\t"
                + result.getNuclSeq().length() + "\t"
                + (protSearch ? result.getProtSeq().length() + "\t" : "")
                + result.getSearchDirection() + "\t"
                + threeDig.format(result.getNats()) + "\t"
                + format.format(result.getBits()) + "\t"
                + threeDig.format(result.getTime()  / 1000.0f));
    }
}
