""" 
Kristina Kolibab, Iqra Rana
Professor Tamon
Assignment 3
"""

import picos as pic
import cvxopt as cvx
import networkx as nx

class Graph:

	# Class initializer
	def __init__(self, filename_, l_, v_, f_):
		self.filename = filename_
		self.graph = [[0 for i in range(v_)] for j in range(v_)]
		self.file_lines = l_
		self.vertices = v_
		self.edges = {}
		self.outputfile = f_
		self.org_graph = [i[:] for i in self.graph]
		self.graph_list = [self.file_lines-2]*0
	
	# Add edges to my graph	
	def addEdge(self, u, v, w):
		self.graph[u][v] = w
		self.graph_list.append(str(u) + "," + str(v))

	# Reads file line by line, creates directed graph
	def readFile(self):
		f = open(self.filename)
		cnt = 1
		n = self.file_lines
		for line in f:
			if(cnt == 1 or cnt == n):
				pass
			else:
				words = line.split()
				u,v = int(words[0]),int(words[2])
				w = [l for l in words[3] if l.isdigit()]
				w = int(''.join(w))
				self.addEdge(u, v, w)
			cnt += 1
		self.org_graph = [i[:] for i in self.graph]
	
	# Used primarily for testing
	def printGraph(self):
		for l in self.graph:
			print(l)

	# Helper method
	def getWeight(self, u, v):
		return self.graph[u][v]
	
	# Helper method for FordFulkerson
	def BFS(self, s, t, parent):
		ROW = len(self.graph)
		visited = [False]*(ROW)
		queue = []
		queue.append(s)
		visited[s] = True
		while queue:
			u = queue.pop(0)
			for ind, val in enumerate(self.graph[u]):
				if visited[ind] == False and val > 0:
					queue.append(ind)
					visited[ind] = True
					parent[ind] = u
		return True if visited[t] else False

	# Method that returns the max flow of a residual graph
	def FordFulkerson(self, source, sink):
		ROW = len(self.graph)
		parent = [-1]*(ROW)
		max_flow = 0
		paths = []
		while self.BFS(source, sink, parent):
			path_flow = float("Inf")
			s = sink
			while(s != source):
				path_flow = min(path_flow, self.graph[parent[s]][s])
				temp = str(parent[s]) + ',' + str(s)	
				paths.append(temp)
				s = parent[s]
			for p in paths:
				if p in self.edges:
					self.edges[p] = self.edges[p] + path_flow
				else:
					self.edges[p] = path_flow
			paths.clear() # don't want duplicates
			max_flow += path_flow
			v = sink
			while(v != source):
				u = parent[v]
				self.graph[u][v] -= path_flow
				self.graph[v][u] += path_flow
				v = parent[v]

		f = open(self.outputfile, "w")
		f.write("digraph {\n")
		for key, value in self.edges.items():
			a, b = key.split(",")
			if key in self.graph_list:
				self.graph_list.remove(key)
			f.write("        " + a + " -> " + b + " [label=")
			f.write('"')
			f.write(str(value))
			f.write('"];\n')
		
		for g in self.graph_list: # for edges that were never used
			a, b = g.split(",")
			f.write("        " + a + " -> " + b + " [label=")
			f.write('"')
			f.write(str(self.graph[int(a)][int(b)]))
			f.write('"];\n')			
		f.write("}")
		f.close()
		return max_flow

	def minCut(self):
		# Make Graph
		G = nx.Graph()
		# Add edges to graph
		ROW = len(self.graph)
		for i in range(ROW):
			for j in range(ROW):
				G.add_edge(i,j)
		# Add weights to edges
		for(i, j) in G.edges():
			G[i][j]['weight'] = self.getWeight(i,j)
	
		# Where the magic happens
		N = self.vertices
		c = {}
		for e in sorted(G.edges(data=True)):			
			capacity = self.graph[e[0]][e[1]]		
			e[2]['capacity'] = capacity
			c[(e[0], e[1])] = capacity
		cc = pic.new_param('c', c)
		s, t = 0, self.vertices-1
		mincut = pic.Problem()
		d = {}
		for e in G.edges():
			d[e] = mincut.add_variable('d[{0}]'.format(e),1)
		p = mincut.add_variable('p', N)
		mincut.add_list_of_constraints(
				[d[i,j] > p[i]-p[j]	
				for (i,j) in G.edges()],
				['i', 'j'], 'edges')
		mincut.add_constraint(p[s] == 1)
		mincut.add_constraint(p[t] == 0)
		mincut.add_constraint(p>0)
		mincut.add_list_of_constraints(
				[d[e] > 0 for e in G.edges()],
				[('e',2)],
				'edges')
		mincut.set_objective('min',
			pic.sum([cc[e]*d[e] for e in G.edges()], [('e',2)], 'edges'))
		mincut.solve(verbose=0, solver='glpk')
		cut = [e for e in G.edges() if abs(d[e].value[0]-1) < 1e-6] 
		S  =[n for n in G.nodes() if abs(p[n].value[0]-1) < 1e-6]
		T  =[n for n in G.nodes() if abs(p[n].value[0]  ) < 1e-6]
		return S, T		
		
def main():
	filename = input("Enter file name: ")
	num_vertices = int(input("Enter the number of vertices: "))
	outputfile = "flow_graph.txt"
	num_lines = sum(1 for line in open(filename))

	g = Graph(filename, num_lines, num_vertices, outputfile)
	g.readFile()

	source = 0; sink = num_vertices-1
	print("Max Flow: %d " % g.FordFulkerson(source, sink))

	S, T = g.minCut()
	fout = open("mincut.txt", "w")	
	fout.write("S Cut: ")
	for s in S:
		fout.write(str(s))
		fout.write(' ')
	fout.write('\n')
	fout.write("T Cut: ")
	for t in T:
		fout.write(str(t))
		fout.write(' ')
	fout.close()

if __name__ == "__main__":
	main()
