#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC
#
#   Generates dispersion models as described at
#   http://yetanothermathprogrammingconsultant.blogspot.com/
#   2019/06/maximum-dispersion.html
#

import sys
import random
import gurobipy as gp
from gurobipy import GRB
from utildev import *

#
#   Solve the p disperson sum problem.
#   Input consists of a list of n points of dimension d as a list of
#   lists.  Compute the distances dij between each point, then formulate
#   the nonconvex MIQP
#   max sum dij*xi*xj
#   s.t.
#   sum xi = k                // k = cardinality of selected points
#   xi binary
#
#   Optionally tighten the model by first linearizing it, then adding
#   a constraint based on the linearization variable (Linearized Model
#   3 in the above URL).
#   Optionally use the objective to use fewer constraints in the linearization
#   process (see the comments for the qlinearize() function below for the
#   details.
#

def pdispersionsum (points, kcard, tighten=False, useobj=False):
    n     = len(points)
    model = gp.Model()
    x     = model.addVars(n, vtype=GRB.BINARY, name="x")
    qobj  = gp.QuadExpr()
#
#   Build the quadratic objective part of the model; just need the upper
#   off diagonals since distances between points are symmetric
#
    for i in range(n):
        for j in range(i+1, n):
            dij = distance(points[i], points[j])
            qobj += dij * x[i] * x[j]
    model.setObjective(qobj, GRB.MAXIMIZE)
#
#   Done with quadratic objective; now add the cardinality constraint.
#
    model.addConstr(gp.quicksum(x) == kcard)
    model.update()
    if tighten:
#
#       Get the quadratic portion of the objective.  Linearize the
#       associated products of binaries.   Then tighten the model
#       by adding the constraint that the sum of the linearization
#       variables must be = kcard*(kcard - 1)/2
#
        zvars = []
        qlinearize(model, useobj, zvars)
        model.addConstr(gp.quicksum(zvars) == kcard * (kcard - 1) / 2, \
                        name="cardinality")
        
    model.update()
    return model
#
#   Takes a model and adds linearization constraints for all
#   products of binaries in the objective.   Optionally returns
#   the list of variables created by the linearization process.
#   Does the standard linearization; no dual arguments to make more
#   compact.   In other words, for xi*xj, create z_ij with
#   zij <= xi                                          (1)
#   zij <= xj                                          (2)
#   zij >= xi + xj - 1                                 (3)
#
#   Since we don't know where the quadratic expression is coming from,
#   is may have separate xi*xj and xj*xi terms.   There is no need
#   to linearize both since zij = zji.   So we avoid creating any
#   duplicates.
#
#   Also, need to adjust the objective.  If original objective is
#   c'x + x'Qx, remove the x'Qx term, keep the c'x term, and add
#   terms for the z variables, where zij takes the same coefficient
#   as xi*xj in the original objective.
#
#   If useobj = True, we use dual arguments based on the objective to
#   create fewer linearization constraints.   If the objective is pushing
#   the z variable down, we only add linearization constraint (1), while
#   we instead add constraints (2) and (3) if the z variable is pushing up.
#
#   TODO: extend to linearize quadratic constraints with products of
#   binaries.
#
def qlinearize(model, useobj, zvars = None):
    qobj     = model.getObjective()
    t        = model.ModelSense
    linpart  = qobj.getLinExpr()
    model.setObjective(linpart)
    dupdict  = {}
    for k in range(qobj.size()):
        qcoef    = qobj.getCoeff(k)
        xi       = qobj.getVar1(k)
        xj       = qobj.getVar2(k)
        key1     = (xi, xj)
        key2     = (xj, xi)
        if key1 in dupdict or key2 in dupdict:
            continue         # this product already linearized; don't duplicate
        else:
            dupdict[key1] = 1   # first time this xi, xj pair encountered.
            dupdict[key2] = 1   # proceed with linearization.            
        qcoeff   = qobj.getCoeff(k)
        zname    = "z_" + xi.VarName + "_" + xj.VarName
        zij      = model.addVar(obj=qcoeff, vtype=GRB.BINARY, name = zname)
        suffix   = xi.VarName + "_" + xj.VarName
        skip1    = False
        skip2    = False
        if useobj:
            down = qcoef*t > 0.0
            if down:
                skip1 = True
            else:               # qobj only contains nonzero Q elements
                skip2 = True
        
        if not skip1:
            cname = "lin1_"+ suffix
            model.addConstr(zij - xi <= 0, name = cname)
            cname = "lin2_"+ suffix
            model.addConstr(zij - xj <= 0, name = cname)
        if not skip2:
            cname = "lin3_"+ suffix
            model.addConstr(xi + xj - zij <= 1, name = cname)
        if zvars != None:
            zvars.append(zij)


#
#   Solve a MIQCP variant of the maximum disperson problem described
#   at the above web site.  Rather than maximizing total distance between
#   points, we maximize the minimum distance, resulting in a set of points
#   more evenly spread out.  As before, input consists of a list of n
#   points of dimension d as a list of lists.  Compute the distances dij
#   between each point, then formulate the nonconvex MIQCP
#
#   max delta
#   s.t.
#   delta <= dij + M(1 - xi*xj)            i < j
#   sum xi = k                // k = cardinality of selected points
#   xi binary
#
#   M can be set to max dij, so as long as we scale our points to have
#   reasonable distances, we won't have any issues with large big M
#   values..
#
#   Optionally tighten the model by replacing
#   
#   delta <= dij + M(1 - xi*xj)            i < j
#   and its associated linearization with the more compact
#
#   delta <= dij + M(1 - xi) + M(1 - xj)
#
            
def betterdispersion (points, kcard, tighten=False):
    n     = len(points)
    model = gp.Model()
    x     = model.addVars(n, vtype=GRB.BINARY, name="x")
    delta = model.addVar(vtype=GRB.CONTINUOUS, name="delta")
#
#   Build the quadratic constraints just need the products xi*xj
#   for i < j since distances between points are symmetric.
#
    model.setObjective(delta, GRB.MAXIMIZE)
    M = 0.0
    distdict = {}
    for i in range(n):
        for j in range(i+1, n):
            dij = distance(points[i], points[j])
            if dij > M:
                M = dij
            distdict[(i,j)] = dij
    if tighten:
        for i in range(n):
            for j in range(i+1, n):
                cname = "dist_" + str(i) + "_" + str(j)
                model.addConstr(delta <= distdict[(i,j)] + 2*M \
                                         - M * x[i] - M * x[j], cname)  
    else:
        for i in range(n):
            for j in range(i+1, n):
                cname = "dist_" + str(i) + "_" + str(j)
                model.addQConstr(delta <= distdict[(i,j)] + M \
                                 - M * x[i] * x[j], cname)
#
#   Done with quadratic constraints; now add the cardinality constraint.
#
    model.addConstr(gp.quicksum(x) == kcard, "cardinality")
    model.update()
    return model
            
        
#
#   Read in points from a file.   The first line is a header line consisting
#   of the number of points to input followed by their dimension
#   The subsequent lines consist of the coordinates of the individual points.
#   Example file of 4 points with dimension 3
#   4   3
#   3.14 2.718 1.6178
#   3.0  2.0   1.0
#   7.14 7.55  7.62
#   19.19 19.59 20.05
#
#   Returns the points as a list of lists
#
def readpointsfromfile(filepath):
    points = []
    if (filepath):
        with open(filepath) as file:
            for i, line in enumerate(file):
                tokens = line.split()
                if (i == 0):             # header line
                    npts = int(tokens[0])
                    dim  = int(tokens[1])
                    if len(tokens) > 2:
                        print("Warning: more than 2 entries in file header: ")
                        print(line)
                        print("Ignoring all but first two entries.")
                    continue
                if len(tokens) != dim:
                    print("Warning: Line ", i, "has ", len(tokens), "entries, ", 
                          "should have ", dim)
                else:       # valid line with coordinates of one point
                    pt = []
                    for j in range(dim):
                        pt.append(float(tokens[j]))
                    points.append(pt)
            
    return points


#
#   Write points from an array file. The first line is a header line consisting
#   of the number of points to input followed by their dimension
#   The subsequent lines consist of the coordinates of the individual points.
#   Example file of 4 points with dimension 3
#   4   3
#   3.14 2.718 1.6178
#   3.0  2.0   1.0
#   7.14 7.55  7.62
#   19.19 19.59 20.05
#
def writepointstofile(filename, points):
    numpts  = len(points)
    dim     = len(points[0])
    ptsfile = open(filename, "w")
    header  = str(numpts) + " " + str(dim)
    ptsfile.write(header + "\n")
    for i in range(numpts):
        line = ""
        for j in range(dim):
            line += str(points[i][j]) + " "
        ptsfile.write(line + "\n")
            
    ptsfile.close()
#
#   Generates num points of dimension dim.   Default domain
#   is [0.0, 1.0] if no domain specified.  Returns points as
#   a list of lists
#
def getrandpoints(num, dim, domain = [0.0, 1.0]):
    points = []
    for i in range(num):
        pt = []
        for j in range(dim):
            pt.append(random.uniform(domain[0], domain[1]))
        points.append(pt)

    return points
