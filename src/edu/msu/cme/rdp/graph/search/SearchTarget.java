package edu.msu.cme.rdp.graph.search;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import edu.msu.cme.rdp.graph.filter.BloomFilter;
import edu.msu.cme.rdp.readseq.SequenceType;

/**
 *
 * @author fishjord
 */
public class SearchTarget {

    private String kmer;
    private int frame;
    private int startState;
    private ProfileHMM forwardHmm;
    private ProfileHMM reverseHmm;
    private BloomFilter filter;
    private boolean prot;

    public SearchTarget(String kmer, int frame, int startState, ProfileHMM forwardHmm, ProfileHMM reverseHmm, BloomFilter filter) {
        this.kmer = kmer;
        this.frame = frame;
        this.startState = startState;
        this.forwardHmm = forwardHmm;
        this.reverseHmm = reverseHmm;
        this.filter = filter;
        this.prot = forwardHmm.getAlphabet() == SequenceType.Protein;
    }

    public int getFrame() {
        return frame;
    }

    public ProfileHMM getForwardHmm() {
        return forwardHmm;
    }

    public ProfileHMM getReverseHmm() {
        return reverseHmm;
    }

    public String getKmer() {
        return kmer;
    }

    public int getStartState() {
        return startState;
    }

    public BloomFilter getFilter() {
        return filter;
    }

    public boolean isProt() {
        return prot;
    }
}
