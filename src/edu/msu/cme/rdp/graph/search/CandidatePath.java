/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.msu.cme.rdp.graph.search;

import edu.msu.cme.rdp.kmer.Kmer;
import java.io.Serializable;
import java.util.*;

/**
 *
 * @author fishjord
 */
public class CandidatePath implements Serializable, Comparable<CandidatePath> {

    List<AStarNode> path = new ArrayList();
    public double score;
    public long generationTime;
    private int iscore;
    int i;
    int k;
    static final long serialVersionUID = 5164991896723207232L;

    private CandidatePath() {
    }

    public CandidatePath(AStarNode goal, Set<Kmer> seenKmers) {
        //score = (int) (goal.score * AStar.INT_SCALE);
        score = goal.score;
        iscore = (int) (goal.score * HMMGraphSearch.INT_SCALE);
        while (goal != null) {
            path.add(goal);
            seenKmers.add(goal.kmer);

            goal = goal.discoveredFrom;
        }

        Collections.reverse(path);
    }

    public CandidatePath subpath(int to) {
        CandidatePath ret = new CandidatePath();
        ret.path = path.subList(0, to);
        //ret.score = (int)(path.get(to).score * AStar.INT_SCALE);
        ret.score = path.get(to).score;
        ret.iscore = (int) (ret.score * HMMGraphSearch.INT_SCALE);

        return ret;
    }

    public AStarNode get(int i) {
        return path.get(i);
    }

    public boolean sameBase(CandidatePath comp, int to) {
        if (comp.path.size() <= to || path.size() <= to) {
            return false;
        }
        for (int index = 0; index <= to; index++) {
            if (!path.get(index).equals(comp.path.get(index))) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String toString() {
        return score + " " + path;
    }

    public void disallowTransition(Map<AStarNode, Set<AStarNode>> disallowedTransitions, int trans) {
        if (path.size() <= trans) {
            return;
        }

        AStarNode node = path.get(trans);

        if (node.discoveredFrom == null) {
            return;
        }

        if (!disallowedTransitions.containsKey(node)) {
            disallowedTransitions.put(node, new HashSet());
        }

        disallowedTransitions.get(node).add(node.discoveredFrom);
    }

    public int length() {
        return path.size();
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (getClass() != obj.getClass()) {
            return false;
        }
        final CandidatePath other = (CandidatePath) obj;
        if (this.path != other.path && (this.path == null || !this.path.equals(other.path))) {
            return false;
        }
        return true;
    }

    @Override
    public int hashCode() {
        int hash = 3;
        hash = 89 * hash + (this.path != null ? this.path.hashCode() : 0);
        return hash;
    }

    public int compareTo(CandidatePath o) {
        int ret = o.iscore - iscore; //Flip because we want the largest first
        if (ret == 0) {
            ret = o.length() - length();
        }

        if (ret == 0) {
            for (int index = 0; index < o.length(); index++) {
                AStarNode n1 = get(index);
                AStarNode n2 = o.get(index);

                if (!n1.equals(n2)) {
                    ret = n1.compareTo(n2);
                    break;
                }
            }
        }

        if (ret == 0) {
            System.err.println("Both paths are the same score length and nodes...that is truly bizarre");
        }

        return ret;
    }
}
