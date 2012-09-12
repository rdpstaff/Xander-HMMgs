/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;
import org.junit.*;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class NextCodonTest {
    /**
     * Test of getCodon method, of class NextCodon.
     */
    @Test
    public void testGetCodon() {

        char[] in = {'a', 'c', 'g', 't'};

        Kmer kmer = new Kmer(in);
        assertEquals(27, kmer.getPart(0));
        assertEquals("acgt", kmer.toString());

        NextCodon codon = new NextCodon(true, new byte[]{0, 1, 2});
        assertEquals('t', codon.getAminoAcid());
        assertEquals(6, codon.getCodon());
        assertEquals("acg", new Kmer(codon.getCodon(), 3).toString());
    }

    @Test
    public void testSanity() {
        int a = 0, c = 1, g = 2, t = 3;
        assertEquals((Character)'s', NextCodon.bacteriaCodonMapping[a][g][t].getAminoAcid());
        assertEquals((Character)'f', NextCodon.bacteriaCodonMapping[t][t][t].getAminoAcid());
        assertEquals((Character)'l', NextCodon.bacteriaCodonMapping[t][t][a].getAminoAcid());
        assertEquals((Character)'g', NextCodon.bacteriaCodonMapping[g][g][t].getAminoAcid());
    }
}
