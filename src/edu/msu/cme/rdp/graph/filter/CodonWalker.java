/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;

/**
 *
 * @author fishjord
 */
public interface CodonWalker {

    public void jumpTo(char[] s);
    public void jumpTo(Kmer kmer, long fwdHash, long rcHash);

    public NextCodon getNextCodon();
    /**
     * @return alternate amino acid for current codon position.
     * If fails, removes current codon (previous codon becomes current).
     * If it backs up to the starting k-mer, it leaves it in that state
     * and returns 0
     *
     */
    public NextCodon getSibCodon();
    public boolean hasMoreCodons();

    public Byte getNextNucl();
    public Byte getSibNucl();
    public boolean hasMoreNucl();

    public long getFwdHash();
    public long getRcHash();

    /**
     *
     * @return the path starting from the char right after the original kmer,
     *  all the way to the current char pointed by pathPtr
     */
    public String getPathString();
}

