require(igraph)

resample <- function(x, ...) {
    x[sample.int(length(x), ...)]
}

ws.graph <- function(n, nei, p) {
    stopifnot(nei < n)
    edge.list <- vector("list", n)
    for (v in 0:(n-1)) {
        edge.end <- union((v + 1:nei) %% n,
                          (v + (-1:-nei)) %% n)
        rewire <- (runif(length(edge.end)) < p)
        edge.end <- edge.end[!rewire]
        rewired <- resample(setdiff(0 : (n-1),
             c(edge.end, v)), sum(rewire))
        edges <- rep(v, 4 * nei)
        edges[c(F, T)] <- c(edge.end, rewired)
        edge.list[[v + 1]] <- edges
    }
    graph(unlist(edge.list))
}