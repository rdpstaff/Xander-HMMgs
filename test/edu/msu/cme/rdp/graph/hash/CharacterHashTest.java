/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package edu.msu.cme.rdp.graph.hash;

import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author wangqion
 */
public class CharacterHashTest {

    public CharacterHashTest() {
    }

    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    /**
     * Test of getHashvalue method, of class CharacterHash.
     */
    @Test
    public void testGetHashvalue() {
        System.out.println("getHashvalue");

        NucleotideHash instance = new NucleotideHash();
        long expResult = 0L;
        long result = instance.getHashvalue(0);
        long lower_result = instance.getHashvalue(0);

        System.err.println("A=" + result + " "+ Long.toBinaryString(result));

        //assertEquals(expResult, result);
        assertEquals(lower_result, result);
        result = instance.getHashvalue(3);
        lower_result = instance.getHashvalue(3);
        System.err.println("T=" + result + " " + Long.toBinaryString(result));
        assertEquals(lower_result, result);

        result = instance.getHashvalue(3);
        lower_result = instance.getHashvalue(3);
        assertEquals(lower_result, result);
    }

    /**
     * Test of getRCHashvalue method, of class CharacterHash.
     */
    @Test
    public void testGetRCHashvalue() {
        /*
        System.out.println("getRCHashvalue");
        CharacterHash instance = new CharacterHash();

        long expResult = 0L;
        long result = instance.getRCHashvalue('A');
        long rev_result = instance.getHashvalue('t');
        assertEquals(rev_result, result);

        result = instance.getRCHashvalue('g');
        rev_result = instance.getHashvalue('c');
        assertEquals(rev_result, result);

        result = instance.getRCHashvalue('u');
        rev_result = instance.getHashvalue('A');
        assertEquals(rev_result, result);

        int x = (int)Math.pow(2, 29);
        System.err.println(x + " " + Integer.toBinaryString(x));
        int wordsize = 32;
        System.err.println("x <<1 =" + (x << 1));
        System.err.println("(x >>> (wordsize-1) =" + (x >>> (wordsize-1) ));

        System.err.println("combined =" + ((x << 1 ) | (x >>> (wordsize-1)) ));

        int m = 17;
        int n = 9;
        System.err.println(m + " " + Integer.toBinaryString(m));
        System.err.println(n + " " + Integer.toBinaryString(n));
        System.err.println("x^ y =" + (m ^ n));
        System.err.println("x^ y ^ y =" + ((m ^ n) ^n));
        System.err.println("x | y =" + (m | n));
        System.err.println("x & y =" + (m & n));

        System.err.println("x >> 1 =" + (-1 >> 1));
*/
    }

}