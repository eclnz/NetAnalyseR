#include <Rcpp.h>
#include <vector>
#include <queue>
#include <limits>

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

    // Implementing Dijkstra's algorithm for each vertex
    for (int src = 0; src < V; src++) {
        std::vector<double> dist(V, std::numeric_limits<double>::max());
        std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;

        dist[src] = 0.0;
        pq.push({0.0, src});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

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

