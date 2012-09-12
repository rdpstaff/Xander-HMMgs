/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.alignment.hmm.ProfileHMM;
import static edu.msu.cme.rdp.alignment.hmm.TSC.*;

/**
 *
 * @author fishjord
 */
public class TestHCostHMM {

    private ProfileHMM hmm;
    private static final char[] mapping = new char[]{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
    private final double[][] mostProbFromState;

    public TestHCostHMM(ProfileHMM hmm) {
        this.hmm = hmm;
        mostProbFromState = new double[128][hmm.M() + 1];

        for (int index = 0; index <= hmm.M(); index++) {
            mostProbFromState['m'][index] = computeCostInternal('m', index);
            mostProbFromState['i'][index] = computeCostInternal('i', index);
            mostProbFromState['d'][index] = computeCostInternal('d', index);
        }
    }

    public double computeHeuristicCost(AStarNode node) {
        return mostProbFromState[node.state][node.stateNo];
    }

    private double computeCostInternal(char prevState, int stateNo) {

        double h = 0;
        double matchTrans = Double.NEGATIVE_INFINITY;
        double insTrans = Double.NEGATIVE_INFINITY;
        double delTrans = Double.NEGATIVE_INFINITY;
        double bestMatchProb = Double.NEGATIVE_INFINITY;
        double bestInsProb = Double.NEGATIVE_INFINITY;

        for (int state = stateNo + 1; state <= hmm.M(); state++) {
            //System.out.println(state + " " + layerScale);

            switch (prevState) {
                case 'm':
                    matchTrans = hmm.tsc(state - 1, MM);
                    insTrans = hmm.tsc(state - 1, MI);
                    delTrans = hmm.tsc(state - 1, MD);
                    break;
                case 'd':
                    matchTrans = hmm.tsc(state - 1, DM);
                    insTrans = Double.NEGATIVE_INFINITY;
                    delTrans = hmm.tsc(state - 1, DD);
                    break;
                case 'i':
                    matchTrans = hmm.tsc(state - 1, IM);
                    insTrans = hmm.tsc(state - 1, II);
                    delTrans = Double.NEGATIVE_INFINITY;
                    break;
                default:
                    throw new RuntimeException("I hate you.");
            }

            bestMatchProb = Float.NEGATIVE_INFINITY;
            bestInsProb = Float.NEGATIVE_INFINITY;

            for (int k = 0; k < hmm.K(); k++) {
                if (hmm.msc(state, k) > bestMatchProb) {
                    bestMatchProb = hmm.msc(state, k);
                }

                if (hmm.isc(state, k) > bestInsProb) {
                    bestInsProb = hmm.isc(state, k);
                }
            }

            matchTrans += bestMatchProb - hmm.getMaxMatchEmission(state);
            delTrans -= hmm.getMaxMatchEmission(state);
            insTrans += bestInsProb;
	    insTrans = Float.NEGATIVE_INFINITY;

            if (insTrans > matchTrans && insTrans > delTrans) {
                h += insTrans;
                prevState = 'i';
                state--;
            } else if (delTrans > matchTrans && delTrans > insTrans) {
                h += delTrans;
                prevState = 'd';
            } else {
                h += matchTrans;
                prevState = 'm';
            }
        }

        return h;
    }
}
