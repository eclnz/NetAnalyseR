#include <Rcpp.h>
#include <vector>
#include <queue>
#include <limits>
#include <algorithm>
#include <utility>

using namespace Rcpp;

// Helper structure to manage edges efficiently
struct Edge {
    int to;
    double weight;
};

typedef std::vector<std::vector<Edge>> Graph;

// Convert adjacency matrix to adjacency list and calculate all-pairs shortest paths
// [[Rcpp::export]]
NumericMatrix dijkstraAllPairs(const NumericMatrix& matrix) {
    int V = matrix.nrow();
    Graph graph(V);

    // Convert matrix to graph (adjacency list)
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            if (matrix(i, j) != 0) {
                graph[i].push_back(Edge{j, matrix(i, j)});
            }
        }
    }

    // Prepare the distance matrix to return
    NumericMatrix distMatrix(V, V);

    // Declared once for the whole all-pairs run rather than once per source, so dist's
    // buffer is allocated a single time and reused.
    typedef std::pair<double, int> QueueEntry;
    std::vector<double> dist(V);
    std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>> pq;

    // Implementing Dijkstra's algorithm for each vertex
    for (int src = 0; src < V; src++) {
        std::fill(dist.begin(), dist.end(), std::numeric_limits<double>::max());
        // std::priority_queue has no clear(), so swap in an empty one to reset it.
        // (This releases its buffer, so only dist actually reuses its storage.)
        std::priority_queue<QueueEntry, std::vector<QueueEntry>, std::greater<QueueEntry>>().swap(pq);

        dist[src] = 0.0;
        pq.push({0.0, src});

        while (!pq.empty()) {
            double d = pq.top().first;
            int u = pq.top().second;
            pq.pop();

            // Lazy deletion: an entry pushed before a shorter path to u was found is
            // stale by the time it is popped, and re-relaxing u's edges would be
            // wasted work.
            if (d > dist[u]) continue;

            for (const auto& edge : graph[u]) {
                int v = edge.to;
                double weight = edge.weight;
                if (dist[u] + weight < dist[v]) {
                    dist[v] = dist[u] + weight;
                    pq.push({dist[v], v});
                }
            }
        }

        // Fill the distance matrix for this source
        for (int i = 0; i < V; i++) {
            distMatrix(src, i) = dist[i];
        }
    }

    return distMatrix;
}

