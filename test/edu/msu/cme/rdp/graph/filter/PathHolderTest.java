/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;
import java.util.Random;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author fishjord
 */
public class PathHolderTest {

    @Test
    public void testPathHolder() {

        Random r = new Random();

        for (int i = 0; i < 100; i++) {
            StringBuilder expected = new StringBuilder();
            int l = r.nextInt(200) + 150;

            PathHolder holder = new PathHolder();
            for (int index = 0; index < l; index++) {
                byte add = (byte) r.nextInt(4);
                holder.push(add);
                expected.append(Kmer.intToChar[add]);
            }

            assertEquals(expected.toString(), new String(holder.toCharArray()));
            StringBuilder buf = new StringBuilder();
            for (int index = 0; index < holder.size(); index++) {
                buf.append(Kmer.intToChar[holder.get(index)]);
            }
            assertEquals(expected.toString(), buf.toString());
        }
    }

    @Test
    public void specificTest() {
        PathHolder path = new PathHolder();
        char[] expected = "aacgttugta".toCharArray();
        String expectedString = new String(expected).replace('u', 't');
        String expectedReverse = new StringBuilder(expectedString).reverse().toString();

        path.init(new Kmer(expectedReverse.toCharArray()));
        assertEquals(expectedReverse, new String(path.toCharArray()));
        assertEquals(expectedReverse, pathToStringViaGet(path));

        path.init(new Kmer(expected));
        assertEquals(expectedString, new String(path.toCharArray()));
        assertEquals(expectedString, pathToStringViaGet(path));

        int substr = expected.length - 1;
        while(path.size() > 0) {
            char exp = expectedString.charAt(substr);
            String e = expectedString.substring(0, substr--);
            char removed = Kmer.intToChar[path.remove()];

            assertEquals(exp, removed);
            assertEquals(e, new String(path.toCharArray()));
            assertEquals(e, pathToStringViaGet(path));
        }
    }

    @Test
    public void pathFromKmerTest() {
        Kmer kmer;
        PathHolder path = new PathHolder();

        kmer = Kmer.randomKmer(30);
        path.init(kmer);
        assertEquals(kmer.toString(), new String(path.toCharArray()));

        kmer = Kmer.randomKmer(60);
        path.init(kmer);
        assertEquals(kmer.toString(), new String(path.toCharArray()));
    }

    private String pathToStringViaGet(PathHolder path) {
        StringBuilder ret = new StringBuilder();
        for(int index = 0;index < path.size();index++) {
            ret.append(Kmer.intToChar[path.get(index)]);
        }

        return ret.toString();
    }
}
