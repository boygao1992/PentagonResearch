########################################
## Generate Reflection Group
########################################

# define a modulo operator for non-integer
def modulo(number,modulus):
    return number - modulus*floor(number / (modulus));

# define a class for reflections
# need to define equivalence of objects for obj in list to work correctly
class Reflection(object):
    def __init__(self, angle):
        self.angle = modulo(angle,pi);

    def returnType(self):
        return self.__class__.__name__;

    def __repr__(self):
        return 'Reflection(' + str(self.angle) + ')'

    def __latex__(self):
        return '\operatorname{Ref}(%s)'.latex(self.angle)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return abs(self.angle - other.angle) in [0, pi]
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(self, other)


# define a class for rotations
# need to define equivalence of objects for obj in list to work correctly
class Rotation(object):
    def __init__(self, angle):
        self.angle = modulo(angle,2*pi);

    def returnType(self):
        return self.__class__.__name__;

    def __repr__(self):
        return 'Rotation(' + str(self.angle) + ')'

    def __latex__(self):
        return '\operatorname{Ref}(%s)'.latex(self.angle)

    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.angle == other.angle
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(self, other)


# group operation
def compose(operation1, operation2):
    pairOfTypes = (operation1.__class__.__name__, operation2.__class__.__name__)
    compositionRules = {
      ('Reflection', 'Reflection') : Rotation(2*(operation1.angle - operation2.angle)),
      ('Rotation', 'Rotation') : Rotation(operation1.angle + operation2.angle),
      ('Rotation', 'Reflection') : Reflection(operation2.angle + operation1.angle/2),
      ('Reflection', 'Rotation') : Reflection(operation1.angle - operation2.angle/2)
    }
    return compositionRules[pairOfTypes]

# note:
# Reflection(0) in set([Reflection(0)]) ==> False
# Reflection(0) in [Reflection(0)] ==> True
# set inclusion does not work properly for defined objects

def generateReflectionGroup(setOfGenerators):
    foundElements = [Rotation(0)]
    newElements = [Rotation(0)]

    while newElements != []:
        temp = newElements;
        newElements = [];

        for elem in temp:
            for gen in setOfGenerators:
                composition = compose(gen, elem)
                if composition not in foundElements:
                    foundElements.append(composition)
                    newElements.append(composition)
    return foundElements

def printReflectionGroup(setOfGenerators):
    counter = 1
    for op in generateReflectionGroup(setOfGenerators):
        print counter, op
        counter += 1

########################################
## Generate Edge Graph and Vertex Graph
########################################

import numpy as np
import networkx as nx

def edgeReflection(r,sAngle):
    if r.__class__.__name__ == 'Reflection':
        return (Reflection(modulo(2*r.angle-sAngle,pi)))
    return (Reflection(modulo(r.angle+sAngle,pi)))

def operationEqual(operation1, operation2):
    if (operation1.__class__.__name__ == operation2.__class__.__name__) and (operation1 == operation2):
        return True
    return False

def edgeMatching(node1, node2,setOfEdges):
    if (node1["edge"] != node2["edge"]):
        return False
    if operationEqual(compose(edgeReflection(node1["operation"],setOfEdges[node1["edge"]]),node1["operation"]),node2["operation"]):
        return True
    return False

def generateEdgeGraph(reflectionGroup,setOfEdges):
    G = nx.Graph()
    index = 0
    for r in reflectionGroup:
        for s in np.arange(len(setOfEdges)):
            index+=1
            G.add_node(index,{"operation":r,"edge":s})
    N = G.number_of_nodes()
    for i in np.arange(1,N):
        for j in np.arange(i+1,N+1):
            if edgeMatching(G.node[i],G.node[j],setOfEdges):
                G.add_edge(i,j)
                G.add_edge(j,i)
                break
    return G

def generateVertexGraph(edgeGraph, num_edge):
    G = nx.Graph()
    for i in edgeGraph.nodes():
        G.add_node(i,{'operation':edgeGraph.node[i]['operation'], 'vertex':edgeGraph.node[i]['edge']})
    for edge in edgeGraph.edges():
        G.add_edge(edge[0],edge[1])
        if (edge[0]) % num_edge == 0:
            G.add_edge(edge[0]+1-num_edge, edge[1]+1-num_edge)
        else:
            G.add_edge(edge[0]+1,edge[1]+1)
    return G

########################################
## Determine Genus of Given Polygon
########################################

# Modified Version
def angleToEdge(setOfAngles):
    setOfEdges = []
    for i in range(len(setOfAngles)):
        if (setOfAngles[i]!=pi):
            if len(setOfEdges)!=0:
                setOfEdges.append(setOfEdges[len(setOfEdges)-1] + (pi - setOfAngles[i]))
            else:
                setOfEdges.append(pi - setOfAngles[i])
    return setOfEdges

# Original
#def angleToEdge(setOfAngles):
#    setOfEdges = []
#    setOfEdges.append(pi - setOfAngles[0])
#    for i in range(1,len(setOfAngles)):
#        setOfEdges.append(setOfEdges[i-1] + (pi - setOfAngles[i]))
#    return setOfEdges

def determineGenus(setOfAngles, detail = False):
    setOfEdges = angleToEdge(setOfAngles)
    setOfGenerators = [Reflection(angle) for angle in setOfEdges]
    reflectionGroup = generateReflectionGroup(setOfGenerators)
    edgeGraph = generateEdgeGraph(reflectionGroup,setOfEdges)
    vertexGraph = generateVertexGraph(edgeGraph, len(setOfEdges))
    cycles = nx.cycle_basis(vertexGraph)
    setOfCycleSize = [len(cycle) for cycle in cycles]
    # G: genus, the number of holes that penetrate the solid
    # E: the number of edges
    # F: the number of faces
    # S: the number of shells
    E = edgeGraph.number_of_edges()
    V = vertexGraph.number_of_nodes() - sum([len(cycle)-1 for cycle in cycles])
    F = len(reflectionGroup)
    S = 1
    if detail:
        print '(E=',E,' V=',V,' F=',F,')'
    return (E-V-F)/2+S


########################################
## Distribution of Genus
########################################

# angle set generator
# Parameters
# ----------
# type : Integer
#     The type of pentagon.
# parameterList : List
#     A set of internal angles. The number of angles is equal to the degree of freedom in the given type of pentagon.
# Returns
#--------
# A.solve_right(b) : Sage Vector
#     The set of all the internal angles in the pentagon.
def angleSetGenerator(type, parameterList):
    matrixList = [
        ## type 6
        [[0,1,0,1,0],
        [0,2,0,0,-1],
        [1,0,1,0,1],
        [1,1,1,1,1],
        [0,1,0,0,0]],
        ## type 7
        [[0,1,0,0,2],
        [0,0,2,1,0],
        [2,1,0,1,0],
        [1,1,1,1,1],
        [0,0,0,0,1]],
        ## type 8
        [[0,2,1,0,0],
        [0,0,0,1,2],
        [2,0,1,1,0],
        [1,1,1,1,1],
        [0,1,0,0,0]],
        ## type 9
        [[2,0,1,0,0],
        [0,0,0,1,2],
        [0,2,1,1,0],
        [1,1,1,1,1],
        [0,0,0,0,1]],
        ## type 10
        [[1,0,0,0,0],
        [0,1,0,0,1],
        [0,1,2,0,0],
        [0,0,1,1,0],
        [0,0,0,0,1]],
        ## type 11
        [[1,0,0,0,0],
        [0,1,1/2,0,0],
        [0,0,1,0,0],
        [0,1,1,1,1],
        [0,0,1,0,1]],
        ## type 12
        [[1,0,0,0,0],
        [0,1,1/2,0,0],
        [0,0,1,0,0],
        [0,1,1,1,1],
        [0,0,1,0,1]],
        ## type 13
        [[0,1,0,0,0],
        [0,0,0,0,1],
        [2,0,0,1,0],
        [0,0,2,1,0],
        [1,0,0,0,0]]
    ]

    vectorList = [
        vector(SR,[pi,0,2*pi,3*pi,-1]), ## type 6, B in (pi/6, pi/2), Singular
        vector(SR,[2*pi,2*pi,2*pi,3*pi,-1]), ## type 7, E in (pi/3, pi), Singular
        vector(SR,[2*pi,2*pi,2*pi,3*pi,-1]), ## type 8, B in (pi/3, pi), Singular
        vector(SR,[2*pi,2*pi,2*pi,3*pi,-1]), ## type 9, E in (3*pi/4,pi), Singular
        vector(SR,[pi/2,pi,2*pi,3*pi/2,-1]), ## type 10, E in (3*pi/10, 7*pi/10)
        vector(SR,[pi/2,pi,-1,5*pi/2,pi]), ## type 11, C in ()
        vector(SR,[pi/2,pi,-1,5*pi/2,pi]), ## type 12, C in ()
        vector(SR,[pi/2,pi/2,2*pi,2*pi,-1]) ## type 13, A in (pi/4,3*pi/4)
    ]

    A = matrix(SR,matrixList[type-6])
    b = vectorList[type-6]
    i,j = 0,0
    while (i<=len(parameterList)-1)&(j<=len(b)-1):
        if b[j] == -1:
            b[j] = parameterList[i]
            i+=1
        j+=1
    return A.solve_right(b)

# farey sequance generator
## Generate all the rational numbers between 0 and 1 with denominators smaller than N
def farey( n, asc=True ):
    ## order: ascending or descending
    if asc:
        a, b, c, d = 0, 1,  1 , n
    else:
        a, b, c, d = 1, 1, n-1, n
    list = []
    list.append(QQ(a/b))
    while (asc and c <= n) or (not asc and a > 0):
        k = int((n + b)/d)
        a, b, c, d = c, d, k*c - a, k*d - b
        list.append(QQ(a/b))
    return list

########################################
## Billard Generator
########################################

def vectorReflection(v, r):
  return 2*(v*r)/(r*r)*r - v;

# generateBilliards
#
# Description
#   Generate billard trajectory.
# Parameters
#   v  : list of vectors, vertices
#   s0 : vector, starting location
#   u0 : vector, starting direction
#   num_iteration : integer, number of iteration
# Returns
#   setOfLocation: list of vectors, location in each iteration
def generateBilliards(v,s0,u0,edge0,num_iteration):
  setOfLocation = [s0]
  s = s0
  u = u0
  edge_index = edge0
  num_edge = len(v)
  terminal = False
  tolerance = 1e-14
  for iteration in range(num_iteration):
    t1_min = oo
    edge_indexs = [x for x in range(num_edge) if x != edge_index]
    for i in edge_indexs:
      ##################################
      ## solve linear system by formulas
      ##################################
      denominator = u[1]*(v[modulo(i+1,num_edge)][0]-v[i][0])-u[0]*(v[modulo(i+1,num_edge)][1]-v[i][1])
      # 0. if solution exists
      if abs(denominator - 0) > tolerance:
        t1 = ((v[modulo(i+1,num_edge)][0]-v[i][0])*(v[i][1]-s[1]) + (v[modulo(i+1,num_edge)][1]-v[i][1])*(s[0]-v[i][0]))/denominator
        t2 = (u[1]*(s[0]-v[i][0]) + u[0]*(v[i][1]-s[1]))/denominator
        # 1. if the solution is not s itself
        if i != edge_index:
          # 2. if the ray hits a vertex
          if abs(t2 - 0) < tolerance or abs(t2 - 1) < tolerance:
            terminal = True
            break
          # 3. if the path is valid (radiating outward and hitting a line segment), select the one with shortest length
          if t1 > 0 and 0 < t2 < 1 and t1 < t1_min:
            edge_index = i
            t1_min = t1
      ###################################################
      ## solve linear system by built-in function 'solve'
      ###################################################
      # 1. if the solution is not s itself
      var('x1,x2')
      equation = [s[0]+u[0]*x1 == v[i][0]+(v[modulo(i+1,num_edge)][0]-v[i][0])*x2,
                  s[1]+u[1]*x1 == v[i][1]+(v[modulo(i+1,num_edge)][1]-v[i][1])*x2]
      solution = solve(equation,x1,x2)
      # 0. if solution exists
      if solution!=[]:
        t1 = solution[0][0].rhs()
        t2 = solution[0][1].rhs()
        # 0. if solution exists
        if t1.is_algebraic():
          # 2. if the ray hits a vertex
          if bool(abs(t2 - 0) < tolerance) or bool(abs(t2 - 1) < tolerance):
            terminal = True
            break
          # 3. if the path is valid (radiating outward and hitting a line segment), select the one with shortest length
          if bool(t1 > 0) and bool(0 < t2 < 1) and bool(t1 < t1_min):
            edge_index = i
            t1_min = t1

    s = s + u*t1_min
    print iteration,':', s
    setOfLocation.append(s)
    u = vectorReflection(u,v[modulo(edge_index+1,num_edge)]-v[edge_index])
    if terminal:
      break
  return setOfLocation
