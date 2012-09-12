/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.readcorrection;

import edu.msu.cme.rdp.graph.search.*;
import java.util.Arrays;
import edu.msu.cme.rdp.alignment.pairwise.ScoringMatrix;
import edu.msu.cme.rdp.graph.filter.CodonWalker;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class NodeEnumerator {

    private static final ScoringMatrix sm = ScoringMatrix.getDefaultNuclMatrix();

    /*static {
        try {
            sm = new ScoringMatrix(ScoringMatrix.class.getResourceAsStream("/data/NUC.4.4"), -100, -1);
        } catch(IOException e) {
            throw new RuntimeException(e);
        }
    }*/

    private int bestMatch = 5;
    private AStarNode next;
    private char[] nextKmer;
    private int index;
    private final boolean forward;
    private final char[] seq;
    private Character nextNucl;

    public NodeEnumerator(boolean forward, char[] seq) {
        this.forward = forward;
        this.seq = seq;
    }

    public Set<AStarNode> enumerateNodes(AStarNode curr, CodonWalker walker) {
        Set<AStarNode> ret = new HashSet();

        //walker.jumpTo(curr.kmer);
        //nextNucl = walker.getNextNucl();

        int remaining;
        if (forward) {
            index = curr.stateNo + 1;
            remaining = seq.length - index;
        } else {
            index = curr.stateNo - 1;
            remaining = index;
        }

        if (nextNucl != null) {

            while (nextNucl != null) {
                nextKmer = makeNextKmer(forward, nextNucl, curr);

                /**
                 * ************************************
                 *
                 * MATCH NODE
                 *
                 *************************************
                 */
                //next = new AStarSearchNode(curr, new char[]{nextNucl}, index, 'm', nextKmer);

                next.score = curr.score + sm.score(nextNucl, seq[index]);
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

                ret.add(next);

                /**
                 * ************************************
                 *
                 * Insert NODE
                 *
                 *************************************
                 */
                //next = new AStarSearchNode(curr, new char[]{nextNucl}, index, 'i', nextKmer);

                next.score = curr.score + sm.getGapOpen();
                next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

                ret.add(next);

                //nextNucl = walker.getSibNucl();
            }

            /**
             * ************************************
             *
             * Delete NODE
             *
             *************************************
             */
            //next = new AStarSearchNode(curr, new char[]{'-'}, curr.stateNo, 'd', curr.kmer);

            next.score = curr.score + sm.getGapOpen();
            next.fval = (int) (HMMGraphSearch.INT_SCALE * (next.score + remaining * bestMatch));

            ret.add(next);
        }

        return ret;
    }

    private static char[] makeNextKmer(boolean forward, char b, AStarNode curr) {
        return null;
        /*char[] kmer;

        if (forward) {
            kmer = Arrays.copyOfRange(curr.kmer, 1, curr.kmer.length + 1);
            kmer[kmer.length - 1] = b;
        } else {
            kmer = new char[curr.kmer.length];

            kmer[0] = b;

            for (int i = 1; i < kmer.length; i++) {
                kmer[i] = curr.kmer[i - 1];
            }
        }

        return kmer;*/
    }
}
