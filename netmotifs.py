# netmotifs.py
import networkx, numpy, scipy.misc, time
import sys

timeout=60.0 # in seconds

def read_graph(fname):
	"""
	Read a graph from a text file.

	Expects the format used in http://moreno.ss.uci.edu/data.html, i.e.
	some metadata followed by a line with the text "DATA:" followed by
	a space delimited matrix of (potentially weighted) adjacency.
	Constructs an unweighted directed graph.
	"""
	with open(fname,'r') as gfile:
		graph = networkx.DiGraph()
		line = gfile.readline()
		while line.strip() != "DATA:":
			line = gfile.readline()
		matlines = gfile.readlines()
		for row in range(len(matlines)):
			graph.add_node(row)
			weights = matlines[row].split()
			for col in range(len(weights)):
				w = float(weights[col])
				if w != 0:
					graph.add_edge(row,col)
					# if (col in graph) and (row in graph[col]):
					# 	graph.add_edge(row,col,double_edge=True)
					# 	graph.add_edge(col,row,double_edge=True)
					# else:
					# 	graph.add_edge(row,col,double_edge=False)
		return graph

def write_edgelist(graph,fname):
	with open(fname,'w') as outf:
		for e in graph.edges():
			outf.write("{} {} 1\n".format(e[0],e[1]))

def link_mask(G,A,B,C):
	# 3 node subgraphs, 6 possible links
	# A->B, B->A
	# A->C, C->A
	# B->C, C->B
	ab, ba = B in G[A], A in G[B]
	ac, ca = C in G[A], A in G[C]
	bc, cb = C in G[B], B in G[C]
	return ab, ba, ac, ca, bc, cb
	# return dict(	ab=(B in G[A]), 
	# 				ba=(A in G[B]),
	# 				ac=(C in G[A]),
	# 				ca=(A in G[C]),
	# 				bc=(C in G[B]),
	# 				cb=(B in G[C]))

subgraphs = [
			#ab		ba		ac		ca		bc		cb
			[False,	True,	False,	False,	True,	False], #B->A,B->C
			[False,	True,	False,	False,	False,	True], # B->A, C->B
			[False,	True,	False,	False,	True,	True], # B->A, B->C, C->B
			[False,	True,	False,	True,	False,	False], # B->A, C->A
			[False,	True,	False,	True,	True,	False], # B->A, B->C, C->A
			[False,	True,	False,	True,	True,	True], # B->A, B->C, C->A, C->B
			[True,	False,	False,	False,	True,	True], # A->B, B->C, C->B
			[True,	True,	False,	False,	True,	True], # A->B, B->A, B->C, C->B
			[True,	False,	False,	True,	True,	False], # A->B, B->C, C->A
			[True,	True,	False,	True,	True,	False], # A->B, B->A, B->C, C->A
			[True,	True,	False,	True,	False,	True], # A->B, B->A, C->A, C->B
			[True,	True,	False,	True,	True,	True], # A->B, B->A, B->C, C->A, C->B
			[True,	True,	True,	True,	True,	True], # ALL
			]

def maskeq(m1,m2):
	for idx in range(len(m1)):
		if m1[idx] != m2[idx]:
			return False
	return True

def count_subgraphs_n3(G):
	counts = numpy.array([0.,]*13)
	lasttime = time.time()
	steps = 0
	laststep = 0
	for A in G:
		for B in G:
			if B == A: 
				continue
			for C in G:
				if (C == B) or (C==A):
					continue
				steps = steps+1
				mask = link_mask(G,A,B,C)
				for idx in range(len(counts)):
					if maskeq(mask,subgraphs[idx]):
						counts[idx] += 1.0
						break
				curtime = time.time()
				if curtime-lasttime > timeout:
					print "Step {} of {}, ({} steps per second)".format(steps,scipy.misc.comb(len(G),3),float(steps-laststep)/float(curtime-lasttime))
					lasttime=curtime
					laststep=steps
	return counts
def randomize_linkswap(G,iters=100,temperature=None,initial_graph=None):
	if initial_graph is None:
		initial_graph = G.copy()
	rv = initial_graph
	for it in range(iters):
		X1 = numpy.random.randint(len(rv))
		if len(rv[X1])==0:
			#No outgoing edges
			continue
		Y1 = rv[X1].keys()[numpy.random.randint(len(rv[X1]))]
		X2 = (X1 + numpy.random.randint(1,len(rv))) % len(rv)
		if len(rv[X2])==0:
			#No outgoing edges
			continue
		Y2 = rv[X2].keys()[numpy.random.randint(len(rv[X2]))]
		if (X1 == Y2) or (X2 == Y1):
			#We picked both halves of a double edge, swapping doesn't do anything
			#OR our swap would introduce a self-loop, so skip
			continue
		# print X1,"->",Y1,":",X2,"->",Y2
		if (Y2 in rv[X1]) or (Y1 in rv[X2]):
			# X1->Y2 or X2->Y1 already exists
			continue
		if (X1 in rv[Y2]) or (X2 in rv[Y1]):
			# We'd be creating a new double edge, so skip
			continue
		# if (rv[X1][Y1]['double_edge']) or (rv[X2][Y2]['double_edge']):
		if (X1 in rv[Y1]) or (X2 in rv[Y2]):
			# if not((rv[X1][Y1]['double_edge']) and (rv[X2][Y2]['double_edge'])):
			if not((X1 in rv[Y1]) and (X2 in rv[Y2])):
				# swapping a non-double edge with a double edge
				continue
			else:
				#remove the selected double edges
				rv.remove_edge(X1,Y1)
				rv.remove_edge(Y1,X1)
				rv.remove_edge(X2,Y2)
				rv.remove_edge(Y2,X2)
				#Add X1->Y2, Y2->X1, and X2->Y1, Y1->X2
				rv.add_edge(X1,Y2)
				rv.add_edge(Y2,X1)
				rv.add_edge(X2,Y1)
				rv.add_edge(Y1,X2)
				continue
		else:
			#NOT double edge, and we've already ensured we're not making a new 
			#double edge in the second if->continue check
			#Remove old edges
			rv.remove_edge(X1,Y1)
			rv.remove_edge(X2,Y2)
			#Add X1->Y2, X2->Y1
			rv.add_edge(X1,Y2)
			rv.add_edge(X2,Y1)
			continue
	return rv

def main(netfname):
	G=read_graph(netfname)
	real_subgraph_counts = count_subgraphs_n3(G)
	random_subgraph_counts = list()
	lasttime = time.time()
	lastiter = 0
	maxswaps = int(len(G.edges())/0.67)
	print "Maximum number of swaps:",maxswaps
	for iteration in range(1000):
		R = randomize_linkswap(G,iters=max(1000,maxswaps))
		random_subgraph_counts.append(count_subgraphs_n3(R))
		curtime = time.time()
		if curtime-lasttime > timeout:
			print "Iteration {} of {} ({} iterations per second)".format(iteration,1000,float(iteration-lastiter)/float(curtime-lasttime))
			lasttime = curtime
			lastiter = iteration
	random_subgraph_counts = numpy.row_stack(random_subgraph_counts)
	print "Real graph counts:"
	print real_subgraph_counts
	print "Random graph counts (mean, then std):"
	print random_subgraph_counts.mean(axis=0)
	print random_subgraph_counts.std(axis=0)
	print "Real graph Z-score"
	zscore = (real_subgraph_counts - random_subgraph_counts.mean(axis=0))/random_subgraph_counts.std(axis=0)
	print  zscore
	print "Significant (p-value < 0.05):"
	for val in numpy.argsort(zscore)[::-1]:
		if zscore[val] < 1.96:
			break
		print val
if __name__ == '__main__':
	main(sys.argv[1])
