/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.graph.abundance;

import edu.msu.cme.rdp.kmer.Kmer;
import edu.msu.cme.rdp.kmer.set.KmerSet;
import edu.msu.cme.rdp.kmer.set.NuclKmerGenerator;
import edu.msu.cme.rdp.readseq.readers.SequenceReader;
import edu.msu.cme.rdp.readseq.readers.Sequence;
import edu.msu.cme.rdp.readseq.utils.IUBUtilities;
import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author fishjord
 */
public class ReadKmerMapper {

    private final int k;
    private final KmerSet<Set<String>> kmerSet;
    private int processedSeqs = 0;

    private long[] val;
    private NuclKmerGenerator kmerGen;
    private Kmer nuclKmer;
    private Set<String> readIds;
    
    public ReadKmerMapper(File contigFile, int k) throws IOException {
        this.k = k;

        List<Sequence> contigSeqs = SequenceReader.readFully(contigFile);

        kmerSet = new KmerSet();

        Sequence seq;
        for (int index = 0; index < contigSeqs.size(); index++) {
            seq = contigSeqs.get(index);
            String seqstr = seq.getSeqString();
            addKmers(seqstr, index);
        }

        kmerSet.printStats();
    }

    private void addKmers(String seqString, int contigIndex) {

        kmerGen = new NuclKmerGenerator(seqString, k);

        while (kmerGen.hasNext()) {
            nuclKmer = kmerGen.next();
            val = nuclKmer.getLongKmers();
            Set<String> kmers = kmerSet.get(val);
            if (kmers == null) {
                kmers = new HashSet();
                kmerSet.add(val, kmers);
            }
        }
    }


    private void processRead(Sequence seq) {
        processRead(seq.getSeqName(), seq.getSeqString());
        processRead(seq.getSeqName(), IUBUtilities.reverseComplement(seq.getSeqString()));
        processedSeqs++;
    }

    private void processRead(String name, String seqString) {

        kmerGen = new NuclKmerGenerator(seqString, k);

        while (kmerGen.hasNext()) {
            val = kmerGen.next().getLongKmers();

            readIds = kmerSet.get(val);

            if (readIds == null) {
                continue;
            }

            readIds.add(name);
        }
    }

    public void printResults(PrintStream out) throws IOException {
        Set<long[]> keys = kmerSet.getKeys();

        kmerGen = new NuclKmerGenerator("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa", k);

        //kmer shouldn't be null yet, so we'll use it to decode the longs
        for(long[] key : keys) {
            out.print(nuclKmer.decodeLong(key));
            for(String readid : kmerSet.get(key)) {
                out.print(" " + readid);
            }
            out.println();
        }
    }

    public static void main(String[] args) throws Exception {
        if (args.length != 3 && args.length != 4) {
            System.err.println("USAGE: ReadKmerMapper <nucl_contig_file> <reads_file> <k> [#threads]");
            System.exit(1);
        }

        final File nuclContigs = new File(args[0]);
        final File readsFile = new File(args[1]);
        final int k = Integer.valueOf(args[2]);
        final int maxThreads;
        final int maxTasks = 25000;

        if(k > 31) {
            System.err.println("k > 31, passing off to the long kmer mapper");
            ReadKmerMapperLongKmer.main(args);
            return;
        }

        if (args.length == 4) {
            maxThreads = Integer.valueOf(args[3]);
        } else {
            maxThreads = Runtime.getRuntime().availableProcessors();
        }

        System.err.println("Starting kmer mapping at " + new Date());
        System.err.println("*  Number of threads:       " + maxThreads);
        System.err.println("*  Max outstanding tasks:   " + maxTasks);
        System.err.println("*  Nucleotide contigs file: " + nuclContigs);
        System.err.println("*  Reads file:              " + readsFile);
        System.err.println("*  Kmer length:             " + k);

        long startTime = System.currentTimeMillis();
        final ReadKmerMapper kmerCounter = new ReadKmerMapper(nuclContigs, k);
        System.err.println("Kmer trie built in " + (System.currentTimeMillis() - startTime) + " ms");

        System.out.println();

        final AtomicInteger processed = new AtomicInteger();
        final AtomicInteger outstandingTasks = new AtomicInteger();


        SequenceReader reader = new SequenceReader(readsFile);
        Sequence seq;

        //ExecutorService service = Executors.newFixedThreadPool(maxThreads);

        startTime = System.currentTimeMillis();
        while ((seq = reader.readNextSequence()) != null) {
            kmerCounter.processRead(seq);
            processed.incrementAndGet();

            if ((processed.get()) % 1000000 == 0) {
                System.err.println("Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");
            }

            /*
             * final Sequence threadSeq = seq;
             *
             * Runnable r = new Runnable() {
             *
             * public void run() { //System.err.println("Processing sequence " +
             * threadSeq.getSeqName() + " in thread " +
             * Thread.currentThread().getName());
             * kmerCounter.processSeq(threadSeq);
             * //System.err.println("Processed count " + processed);
             * //System.err.println("Outstanding count count " +
             * outstandingTasks); processed.incrementAndGet();
             * outstandingTasks.decrementAndGet(); } };
             *
             * outstandingTasks.incrementAndGet(); service.submit(r);
             *
             * //System.err.println("Submitting " + threadSeq.getSeqName() + ",
             * outstanding tasks= " + outstandingTasks);
             *
             * while (outstandingTasks.get() >= maxTasks);
             *
             * if ((processed.get() + 1) % 1000000 == 0) {
             * System.err.println("Processed " + processed + " sequences in " +
             * (System.currentTimeMillis() - startTime) + " ms"); }
             */
        }

        reader.close();

        //service.shutdown();
        //service.awaitTermination(1, TimeUnit.DAYS);

        System.err.println("Processed " + processed + " sequences in " + (System.currentTimeMillis() - startTime) + " ms");

        kmerCounter.printResults(System.out);
        //System.err.println("Unique kmers in contigs: " + kmerCounter.trie.uniqueWords());
        System.err.println("Processing complete");
    }
}
