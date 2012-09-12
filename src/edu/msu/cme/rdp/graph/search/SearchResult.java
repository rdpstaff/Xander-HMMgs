package edu.msu.cme.rdp.graph.search;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
/**
 *
 * @author fishjord
 */
public class SearchResult {
    public static enum SearchDirection { left, right };
    
    private String kmer;
    private String nuclSeq;
    private String protSeq;
    private long time;
    private double nats;
    private double bits;
    private SearchDirection searchDirection;
    private int mpos;
    private String alignSeq;

    public SearchResult(String kmer, String nuclSeq, String alignSeq, String protSeq, SearchDirection dir, int mpos, double nats, double bits, long time) {
        this.bits = bits;
        this.nats = nats;
        this.time = time;
        this.mpos = mpos;

        this.kmer = kmer;
        this.nuclSeq = nuclSeq;
        this.protSeq = protSeq;
        this.searchDirection = dir;
        this.alignSeq = alignSeq;
    }

    public String getKmer() {
        return kmer;
    }

    public String getNuclSeq() {
        return nuclSeq;
    }

    public String getProtSeq() {
        return protSeq;
    }

    public long getTime() {
        return time;
    }

    public double getBits() {
        return bits;
    }

    public double getNats() {
        return nats;
    }

    public SearchDirection getSearchDirection() {
        return searchDirection;
    }

    public int getModelPosition() {
        return mpos;
    }

    public String getAlignSeq() {
        return alignSeq;
    }
}
