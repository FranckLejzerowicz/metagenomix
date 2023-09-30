# ----------------------------------------------------------------------------
# Copyright (c) 2022, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from collections import defaultdict


class Graph(object):
    """https://www.geeksforgeeks.org/find-paths-given-source-destination/"""
    def __init__(self, vertices):
        self.V = vertices
        self.graph = defaultdict(list)
        self.paths = defaultdict(list)

    def add_edge(self, u, v):
        self.graph[u].append(v)

    def print_paths_util(self, u, d, visited, path):
        visited[u] = True
        path.append(u)
        if u == d:
            self.paths[path[-1]].append(list(path))
        else:
            for i in self.graph[u]:
                if not visited[i]:
                    self.print_paths_util(i, d, visited, path)
        path.pop()
        visited[u] = False

    def print_paths(self, s, d):
        visited = [False] * self.V
        path = []
        self.print_paths_util(s, d, visited, path)
