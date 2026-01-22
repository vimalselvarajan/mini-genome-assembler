import networkx as nx

def parse_overlap_graph_nx(filename):
    G = nx.DiGraph()
    with open(filename) as f:
        for line in f:
            src, dst, ovl = line.strip().split()
            ovl = int(ovl)
            G.add_edge(src, dst, overlap=ovl)
    return G

def is_1_in_1_out(G, node):
    return G.in_degree(node) == 1 and G.out_degree(node) == 1

def maximal_non_branching_paths_nx(G):
    paths = []
    visited_edges = set()

    # Non 1-in-1-out starts
    for v in G.nodes():
        if not is_1_in_1_out(G, v):
            if G.out_degree(v) > 0:
                for w in G.successors(v):
                    edge = (v, w)
                    if edge in visited_edges:
                        continue
                    path = [v, w]
                    visited_edges.add(edge)
                    # extend path while successor is 1-in-1-out
                    while is_1_in_1_out(G, w):
                        successors = list(G.successors(w))
                        if not successors:
                            break
                        u = successors[0]
                        next_edge = (w, u)
                        if next_edge in visited_edges:
                            break

                        if not is_1_in_1_out(G, u):
                            if G.out_degree(u) == 0:
                                path.append(u)  # Explicitly add dead end as last node
                            break  # Do not add nxt! It starts a new unitig
                            
                        path.append(u)
                        visited_edges.add(next_edge)
                        w = u
                    paths.append(path)

    # Cycles: all nodes 1-in-1-out
    for v in G.nodes():
        if is_1_in_1_out(G, v):
            for w in G.successors(v):
                edge = (v, w)
                if edge in visited_edges:
                    continue
                # Begin cycle from this edge
                cycle = [v, w]
                visited_edges.add(edge)
                current = w
                while True:
                    if not is_1_in_1_out(G, current):
                        break
                    successors = list(G.successors(current))
                    if len(successors) != 1:
                        break
                    nxt = successors[0]
                    next_edge = (current, nxt)
                    if next_edge in visited_edges:
                        break
                    cycle.append(nxt)
                    visited_edges.add(next_edge)
                    current = nxt
                    if current == v:
                        paths.append(cycle)
                        break

    return paths

def write_unitigs_with_overlaps(paths, G, outfilename):
    with open(outfilename, 'w') as f:
        for i, path in enumerate(paths):
            f.write(f"START OF UNITIG {i} {path[0]}\n")
            if len(path) > 1:
                for u, v in zip(path[:-1], path[1:]):
                    ovl = G.edges[u, v]['overlap']
                    f.write(f"{v} {ovl}\n")
            f.write(f"END OF UNITIG {i}\n\n")

if __name__ == "__main__":
    G = parse_overlap_graph_nx("overlaps.txt")
    paths = maximal_non_branching_paths_nx(G)
    write_unitigs_with_overlaps(paths, G, "unitigs.txt")
    print(f"Done! Found and wrote {len(paths)} unitigs to 'unitigs.txt'.")
