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
package edu.msu.cme.rdp.graph;

import edu.msu.cme.rdp.graph.search.AStarNode;
import edu.msu.cme.rdp.kmer.Kmer;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;

/**
 *
 * @author fishjord
 */
public class Graph implements Serializable {

    protected static class GraphNode implements Serializable {

        private final int hash;
        public final Kmer kmer;
        private final int mpos;
        private final char state;

        public GraphNode(Kmer kmer) {
            this.kmer = kmer;
            mpos = 0;
            state = 0;
            hash = hashMe(this);
        }

        public GraphNode(Kmer kmer, int mpos, char state) {
            this.kmer = kmer;
            this.mpos = mpos;
            this.state = state;
            hash = hashMe(this);
        }

        @Override
        public boolean equals(Object obj) {
            if (obj == null) {
                return false;
            }
            if (getClass() != obj.getClass()) {
                return false;
            }
            final GraphNode other = (GraphNode) obj;

            if (this.mpos != other.mpos) {
                return false;
            }

            if (this.state != other.state) {
                return false;
            }

            if (this.kmer != other.kmer) {
                return false;
            }
            return true;
        }

        private static int hashMe(GraphNode node) {
            int hash = 3;
            hash = 89 * hash + node.kmer.hashCode();
            hash = 89 * hash + node.mpos;
            hash = 89 * hash + node.state;
            return hash;
        }

        @Override
        public int hashCode() {
            return hash;
        }
    }
    protected Map<GraphNode, Set<GraphNode>> edges = new HashMap();
    protected Map<GraphNode, GraphNode> nodes = new HashMap(); //Don't ask
    protected boolean combined;
    private final int k;

    public Graph(boolean combined, int k) {
        this.combined = combined;
        this.k = k;
    }

    private List<Set<GraphNode>> getComponents() {
        List<Set<GraphNode>> ret = new ArrayList();
        if (nodes.isEmpty()) {
            return ret;
        }


        /*
         * Replace all the diedges with edges
         */
        Set<GraphNode> allNodes = new HashSet(nodes.keySet());
        Map<GraphNode, Set<GraphNode>> weakEdges = new HashMap();
        for (GraphNode ls : edges.keySet()) {
            if (!weakEdges.containsKey(ls)) {
                weakEdges.put(ls, new HashSet());
            }

            weakEdges.get(ls).addAll(edges.get(ls));

            for (GraphNode rs : edges.get(ls)) {
                if (!weakEdges.containsKey(rs)) {
                    weakEdges.put(rs, new HashSet());
                }

                weakEdges.get(rs).add(ls);
            }
        }

        while (!allNodes.isEmpty()) {
            GraphNode next = allNodes.iterator().next(); //Graph any node, don't care
            Set<GraphNode> comp = new HashSet();

            connectedTo(next, weakEdges, comp);

            ret.add(comp);
            allNodes.removeAll(comp);
        }


        return ret;
    }

    private static void connectedTo(GraphNode node, Map<GraphNode, Set<GraphNode>> edges, Set<GraphNode> comp) {
        if (comp.contains(node)) {
            return;
        }

        comp.add(node);
        for (GraphNode next : edges.get(node)) {
            connectedTo(next, edges, comp);
        }
    }

    public void writeDot(File outStem) throws IOException {
        System.out.println("Graph nodes: " + nodes.size());
        System.out.println("Edges: " + edges.size());
        int compCount = 0;
        for (Set<GraphNode> comp : getComponents()) {
            compCount++;
            System.out.println("Component " + compCount + " has " + comp.size() + " nodes");

            PrintWriter out = new PrintWriter(new File(outStem.getAbsolutePath() + "comp_" + compCount + ".gv"));
            try {
                out.println("digraph component_" + compCount + " {");
                for (GraphNode ls : comp) {
                    String lskmer = ls.kmer.toString();
                    out.println("/* " + lskmer + " */");
                    for (GraphNode rs : edges.get(ls)) {
                        out.println("\t" + lskmer + " -> " + rs.kmer + ";");
                    }
                    out.println();
                }
                out.println("}");

            } finally {
                out.close();
            }
        }
    }

    public void extend(AStarNode node) {
        extendInternal(node);
    }

    private GraphNode extendInternal(AStarNode node) {
        if(node.state == 'd' && !combined) {
            return extendInternal(node.discoveredFrom);
        }

        GraphNode gn;

        if (combined) {
            gn = new GraphNode(node.kmer, node.stateNo, node.state);
        } else {
            gn = new GraphNode(node.kmer);
        }

        if (nodes.containsKey(gn)) {
            gn = nodes.get(gn); //Don't store duplicated information
        } else {
            nodes.put(gn, gn);
            edges.put(gn, new HashSet());
        }

        if (node.discoveredFrom != null) {
            GraphNode back = extendInternal(node.discoveredFrom);
            edges.get(back).add(gn);
        }

        return gn;
    }
}
