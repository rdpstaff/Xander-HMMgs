/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.graph.filter.PathHolder;
import edu.msu.cme.rdp.kmer.Kmer;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class AStarSearchNodeTest {

    public AStarSearchNodeTest() {
    }

    @Test
    public void testComparison() {
        //AStarSearchNode discoveredFrom, long kmer, long fwdHash, long rcHash, int stateNo, char state
        AStarNode node1 = new AStarNode(null, new Kmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'm');
        AStarNode node2 = new AStarNode(null, new Kmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'd');
        AStarNode node3 = new AStarNode(null, new Kmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'i');
        AStarNode node4 = new AStarNode(null, new Kmer(new char[]{'a', 'c', 'g'}), 0, 0, 15, 'm');
        AStarNode node5 = new AStarNode(null, new Kmer(new char[]{'a', 'a', 'a'}), 0, 0, 15, 'm');
        AStarNode node6 = new AStarNode(null, new Kmer(new char[]{'c', 'g', 'a'}), 0, 0, 16, 'm');
        AStarNode node7 = new AStarNode(null, new Kmer(new char[]{'t', 'a', 'c'}), 0, 0, 50, 'm');

        assertTrue(node1.equals(node1));
        assertTrue(node1.equals(node5));
        assertTrue(node5.equals(node1));

        assertFalse(node1.equals(node2));
        assertFalse(node1.equals(node3));
        assertFalse(node1.equals(node4));
        assertFalse(node1.equals(node6));
        assertFalse(node1.equals(node7));

        assertTrue(node1.compareTo(node2) < 0);
        assertTrue(node1.compareTo(node1) == 0);
        assertTrue(node1.compareTo(node7) < 0);
        assertTrue(node7.compareTo(node1) > 0);

        node1.fval = -1;
        node6.fval = -6;

        node2.fval = -10;

        assertTrue(node2.compareTo(node6) > 0);
        assertTrue(node6.compareTo(node2) < 0);

        assertTrue(node1.compareTo(node6) < 0);
        assertTrue(node6.compareTo(node1) > 0);
    }
}
