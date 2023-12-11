# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict


class Graph(object):
    """
    This class represents a directed graph using adjacency list representation
    https://www.geeksforgeeks.org/find-paths-given-source-destination/"""
    def __init__(self, vertices):
        # No. of vertices
        self.V = vertices
        # default dictionaries to store graph and paths through graph
        self.graph, self.paths = defaultdict(list), defaultdict(list)

    def add_edge(self, u, v):
        """function to add an edge to graph"""
        self.graph[u].append(v)

    def print_paths_util(self, u, d, visited, path):
        # Mark the current node as visited and store in path
        visited[u] = True
        path.append(u)
        if u == d:
            # If current vertex is same as destination, collect the current path
            self.paths[path[-1]].append(list(path))
        else:
            # If not, recurse for all the vertices adjacent to this vertex
            for i in self.graph[u]:
                if not visited[i]:
                    self.print_paths_util(i, d, visited, path)
        # Remove current vertex from path and mark it as unvisited
        path.pop()
        visited[u] = False

    def print_paths(self, s, d):
        visited = [False] * self.V
        path = []
        self.print_paths_util(s, d, visited, path)
