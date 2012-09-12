/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph;

import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.trie.KmerGenerator;
import edu.msu.cme.rdp.readseq.readers.SeqReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import java.io.File;
import java.io.IOException;

/**
 *
 * @author fishjord
 */
public class PostProcessContigs {

    public static void main(String[] args) throws IOException {
        if(args.length != 2) {
            System.err.println("USAGE: PostProcessContigs <nucl_contigs> <k>");
            System.exit(1);
        }

        SeqReader reader = new SequenceReader(new File(args[0]));
        int k = Integer.valueOf(args[1]);

        Sequence seq;
        Graph graph = new Graph(false, k);

        while((seq = reader.readNextSequence()) != null) {
            AStarNode node = null;
            for(char[] kmer : KmerGenerator.getKmers(seq.getSeqString(), k)) {
    //public AStarNode(AStarNode discoveredFrom, long kmer, long fwdHash, long rcHash, int stateNo, char state) {
                AStarNode next = new AStarNode(node, new Kmer(kmer), 0, 0, 0, 'm');
                node = next;
            }

            graph.extend(node);
        }

        graph.writeDot(new File("ppg"));

        reader.close();
    }

}
