import java.text.DecimalFormat;
import java.util.*;

public class MST {
    MST() {

    }
    static double getTotalEdgeWeight(Graph g) {
        double weight = 0;
        for (int i = 0; i < g.numVertices; i++) {
            for (int x = 0; x < g.numVertices; x++) {
                if (g.isEdge(i, x)) {
                    weight += g.weight(i,x);
                }
            }
        }
        return weight/2;
    }
    static Graph getRandomGraph (int n) {
        Graph g = new MatrixGraph(n, Graph.UNDIRECTED_GRAPH);
        double[][] coords = new double[n][2];
        double x;
        double y;
        double x2;
        double y2;
        Random rand = new Random();
        DecimalFormat df2 = new DecimalFormat( "#,###,###,##0.00" );
        for (int i = 0; i < g.numVertices(); i++) {
            x = rand.nextFloat();
            y = rand.nextFloat();
            x = Double.parseDouble(df2.format(x));
            y = Double.parseDouble(df2.format(y));
            coords[i][0] = x;
            coords[i][1] = y;

        }
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                if (j != k) {
                    if (!g.isEdge(j,k)) {
                        x = coords[j][0];
                        y = coords[j][1];
                        x2 = coords[k][0];
                        y2 = coords[k][1];
                        double result = x - x2;
                        double result2 = y - y2;
                        result = result * result;
                        result2 = result2 * result2;
                        result = result + result2;
                        double weight = Double.parseDouble(df2.format(Math.sqrt(result)));
                        g.addEdge(j, k, weight);
                    }
                }
            }
        }
        return g;
    }
    static Graph getBaseTree (Graph g) {
        Graph spanningTree = new MatrixGraph(g.numVertices, Graph.UNDIRECTED_GRAPH);
        int startVertex = 0;
        int cur;
        Queue<Integer> queue = new LinkedList<>();
        List<Integer> visited = new ArrayList<>();
        queue.add(startVertex);
        visited.add(startVertex);
        while (!queue.isEmpty()) {
            cur = queue.remove();
            int[] outNeighbours = g.outNeighbours(cur);
            boolean visitedYes = false;
            for (int i = 0; i < outNeighbours.length; i++) {
                for (int j = 0; j < visited.size(); j++) {
                    if (visited.get(j) == outNeighbours[i]) {
                        visitedYes = true;
                        break;
                    }
                }
                if (!visitedYes) {
                    spanningTree.addEdge(cur, outNeighbours[i], g.weight(cur, outNeighbours[i]));
                    visited.add(outNeighbours[i]);
                    queue.add(outNeighbours[i]);
                }
            }
        }
        return spanningTree;
    }

    static class Node implements Comparable<Node> {
        int v;
        double distance;

        public Node(int v, double distance)
        {
            this.v = v;
            this.distance = distance;
        }

        @Override public int compareTo(Node n)
        {
            if (this.distance <= n.distance) {
                return -1;
            }
            else {
                return 1;
            }
        }
    }

    static Edge longestEdgeOnPath (Graph g, int source, int dest) {
        Edge longestEdge = new Edge(0, 0 , 0);
        if (source == dest) {
            return null;
        }
        int cur;
        boolean[] solved = new boolean[g.numVertices];
        double[] distance = new double[g.numVertices];
        int[] visited = new int[g.numVertices];
        Arrays.fill(visited, -1);
        PriorityQueue<Node> queue = new PriorityQueue<>();
        for (int v = 0; v < g.numVertices; v++) {
            solved[v] = false;
            distance[v] = Double.POSITIVE_INFINITY;
        }
        distance[source] = 0;
        visited[source] = 0;
        queue.add(new Node(source, 0));

        while (!queue.isEmpty()) {
            Node n = queue.remove();
            cur = n.v;
            //System.out.println(n.distance + " " + n.v);
            solved[cur] = true;
            int[] Neighbours = g.outNeighbours(cur);
            //System.out.println(Arrays.toString(Neighbours));
            for (int i = 0; i < Neighbours.length; i++) {
                int y = Neighbours[i];
                if (y == dest) {
                    double newDistance = distance[cur] + g.weight(cur, y);
                    distance[y] = newDistance;
                    visited[y] = cur;
                    queue.clear();
                    break;
                } else if (!solved[y]) {
                    double newDistance = distance[cur] + g.weight(cur, y);
                    if (newDistance < distance[y]) {
                        if (queue.contains(new Node(y, distance[y]))) {
                            queue.remove(new Node(y, distance[y]));
                            queue.add(new Node(y, newDistance));
                        } else {
                            queue.add(new Node(y, newDistance));
                        }
                        distance[y] = newDistance;
                        visited[y] = cur;
                    }
                }
            }
        }
        //System.out.println(Arrays.toString(distance));
        //System.out.println(Arrays.toString(visited));
        double length = 0;
        int index = 0;
        for (int x = 0; x < distance.length; x++) {
            if (distance[x] != Double.POSITIVE_INFINITY) {
                if (distance[x] > length) {
                    length = distance[x];
                    index = x;
                }
            }
        }
        longestEdge.x = visited[index];
        longestEdge.y = index;
        longestEdge.w = g.weight(longestEdge.x, longestEdge.y);
        return longestEdge;
    }
    static Graph getMST (Graph g) {
        Graph spanningTree = getBaseTree(g);

        for (int j = 0; j < g.numVertices; j++) {
            int[] Neighbours = g.outNeighbours(j);
            for (int k = 0; k < Neighbours.length; k++) {
                if (spanningTree.isEdge(j,Neighbours[k])) {
                    Edge e = longestEdgeOnPath(spanningTree, j, Neighbours[k]);
                    if (g.weight(j, Neighbours[k]) < e.w) {
                        spanningTree.deleteEdge(e.x, e.y);
                        spanningTree.addEdge(j, Neighbours[k], g.weight(j, Neighbours[k]));
                    }
                }
            }
        }
        return spanningTree;
    }
    public static void main(String[] args) {
        Graph g = GraphOfEssex.getGraph();
        System.out.println(getTotalEdgeWeight(g));
        Graph spanningTree = getMST(g);
        System.out.println(getTotalEdgeWeight(spanningTree));
        double averageEdgeWeight = 0;
        double averageEdgeWeightMST = 0;
        for (int i = 0; i < 100; i++) {
            Graph randomGraph = getRandomGraph(50);
            averageEdgeWeight += getTotalEdgeWeight(randomGraph);
            Graph MST = getMST(randomGraph);
            averageEdgeWeightMST += getTotalEdgeWeight(MST);
        }
        System.out.println(averageEdgeWeight/100);
        System.out.println(averageEdgeWeightMST/100);
    }
}

