/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.cli;

import edu.msu.cme.rdp.graph.filter.BloomFilter;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.io.PrintStream;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author fishjord
 */
public class BloomFilterStats {
    
    public static void printStats(BloomFilter filter, PrintStream out) {
        

	long n = filter.getUniqueKmers();
	long m = (long)Math.pow(2, filter.getHashSizeLog2());
	int k = filter.getHashCount();

        //(1-e^(-k*((n+.5)/(m-1))))^k
        double falsePositiveRate = Math.pow((1 - Math.pow(Math.E, -k * ((n + .5) / (m - 1)))), k);
        
        out.println("Bloom filter created on:       " + filter.getCreatedOn());
        out.println("Serializable id:               " + BloomFilter.serialVersionUID);
        out.println();
        out.println("Bloom filter size log 2:       " + filter.getHashSizeLog2());
        out.println("Bloom filter size (bits) (m):  " + m);
        out.println();
        out.println("Number of bitsets:             " + filter.getNumBitsets());
        out.println("Bitset size (bits):            " + filter.getBitsetSize());
        out.println("Bitset size log2:              " + filter.getBitsetSizeLog2());
        out.println();
        out.println("Number of hashes (k):          " + filter.getHashCount());
        out.println("Hash function name:            " + filter.getHasherClassName());
        out.println();
        out.println("Bitset Mask:                   " + StringUtils.leftPad(Long.toBinaryString(filter.getBitsetMask()), 64, '0'));
        out.println("Hash mask:                     " + StringUtils.leftPad(Long.toBinaryString(filter.getHashMask()), 64, '0'));
        out.println();
        out.println("Kmer length:                   " + filter.getKmerSize());
        out.println();
        out.println("Total kmers in bloom filter:   " + filter.getTotalKmers());
        out.println("Total strings inserted:        " + filter.getTotalStrings());
        out.println("Total unique kmers:            " + filter.getUniqueKmers());
        out.println("Predicted false positive rate: " + falsePositiveRate);
    }
    
    public static void main(String[] args) throws Exception {
        if(args.length != 1) {
            System.err.println("USAGE: BloomFilterStats <bloom_filter>");
            System.exit(1);
        }
        
        File bloomFile = new File(args[0]);
        
        ObjectInputStream ois = ois = new ObjectInputStream(new BufferedInputStream(new FileInputStream(bloomFile)));
        BloomFilter filter = (BloomFilter) ois.readObject();
        ois.close();
        
        printStats(filter, System.out);
    }
}
