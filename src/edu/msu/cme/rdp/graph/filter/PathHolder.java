/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.filter;

import edu.msu.cme.rdp.kmer.Kmer;
import java.util.Arrays;

/**
 *
 * @author fishjord
 */
public class PathHolder {

    long[] path = new long[50];
    private int ptr = 0;
    private int cap = 32;
    private int size = 0;

    public void init(Kmer kmer) {
        int k = kmer.length();
        //System.err.println("k=" + kmer.length() + ", packedlength: " + kmer.packedLength());

        for(int index = kmer.packedLength() - 1;index >= 0 ;index--) {
            path[index] = kmer.getPart(index);
        }
        ptr = (int)Math.floor(k / 32.0);
        cap = 32 - k % 32;
        size = k;
    }

    public void push(byte twobits) {
        if (cap == 0) {
            ptr++;
            cap = 32;

            if (ptr >= path.length) {
                path = Arrays.copyOf(path, ptr + 50);
            }
            path[ptr] &= 0;
        }

        path[ptr] = (path[ptr] << 2) | (twobits & 3);
        cap--;
        size++;
    }

    public byte remove() {
        if (cap == 32) {
            cap = 0;
            ptr--;
        }

        byte ret = (byte)(path[ptr] & 3);
        path[ptr] = path[ptr] >>> 2;
        size--;
        cap++;

        return ret;
    }

    public int size() {
        return size;
    }

    public byte get(int index) {
        if (index >= size) {
            throw new ArrayIndexOutOfBoundsException(index);
        }

        int targetptr = index / 32;
        int targetcap;

        if(targetptr == ptr) {
            targetcap = (31 - cap) - (index % 32);
        } else {
            targetcap = 31 - (index % 32);
        }

        return (byte) ((path[targetptr] >>> (targetcap * 2)) & 3);
    }

    public char[] toCharArray() {
        char[] ret = new char[size];

        int currptr = ptr;
        int currcap = 32 - cap;
        long val = path[ptr];

        for (int index = size - 1; index >= 0; index--) {
            if (currcap == 0) {
                currptr--;
                currcap = 32;
                val = path[currptr];
            }

            ret[index] = Kmer.intToChar[(int) (val & 3)];
            val = val >>> 2;
            currcap--;
        }

        return ret;
    }
}
