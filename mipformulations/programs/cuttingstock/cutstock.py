#
# One dimensional Cutting stock model formulation and solution.
# Can input the data in the first three arguments, or by file specified
# in the fourth argument. Can generate either a compact but weak
# formulation or the stronger pattern based formulation with the fifth argument.
#
import gurobipy as gp
from gurobipy import GRB

import math
import time

PATTERNS = 1
COMPACT  = 2
SUCCEED  = 1
FAIL     = 0

def cutstock(stockwidth=None, cutwidths=None, cutdemands = None, \
             inputfile=None, modeltype=PATTERNS, verbose=False):
    #
    #   Help function info
    #
    """
    Create a cutting stock model and solve it with simple column generation.

    Arguments:
    stockwidth     Stock roll widths
    cutwidths      List of product roll widths
    cutdemands     List of demand for each product roll.  Must align with 
                   cutwidths array, i.e., cutdemands[i] is the demand for 
                   product roll with width cutwidths[i]
    inputfile      Instead of passing the 3 previous arguments in the 
                   function, specify the problem data in a file.  The 
                   first line specifies the stock roll width. Each subsequent
                   line specifies a product roll width and the associated
                   product roll demand.   Data must be specified either 
                   via this argument or the preceding 3 arguments.
    modeltype      Set to PATTERN (the default) to invoke column generation,
                   or COMPACT to generate the weaker but compact MILP model.
    verbose        Set to True to obtain additional column generation algorithm
                   output.
    Returns:       A model if either the final restricted master or the 
                   compact formulation, depending on the modeltype setting."""
    if inputfile != None:
        if stockwidth != None or cutwidths != None or cutdemands != None:
            print("Warning file specified but non empty data provided ")
            print("directly for width or demand data.  Using file data.")
        with open(inputfile) as file:
            #
            # File format
            # Line 1: Stock width
            # Each subsequent line contains a cutwidth length and the
            # associated demand for that cutwidth.
            #
            cutwidths  = []
            cutdemands = []
            for i, line in enumerate(file):
                tokens = line.split()
                if i == 0:
                    if len(tokens) != 1:
                        print("Error: Single stock width not properly ",\
                               "specified in first line.   Exiting.")
                        return None
                    stockwidth = float(tokens[0])
                else:
                    if len(tokens) != 2:
                        print("Error: lines after first one must specify ",\
                              "a cut width and cutdemand pair.  Exiting.")
                        return None
                    cutwidths.append(float(tokens[0]))
                    cutdemands.append(float(tokens[1]))
    #
    # Data is set; do some sanity checking.
    #
    numwidths = len(cutwidths)
    if numwidths != len(cutdemands):
        print("Error: mismatch in cut widths and cut demand array lengths: ", \
              numwidths, " and ", len(cutdemand), ". Exiting.")
        return None
    for j in range(len(cutwidths)):
        if cutwidths[j] <= 0 or cutdemands[j] < 0:
            print("Error: bogus value in width demand pair ", j), \
                  ": (", cutwidths[j], ",", cutdemands[j], ")"
            print("Exiting.")
            return None

    if modeltype == COMPACT:
        return compact_formulation(stockwidth, cutwidths, cutdemands, verbose)
        

        
    starttime = time.time()
    rm = init_restricted_master(stockwidth, cutwidths, cutdemands)
    if rm == None:
        print("Unable to create initial restricted master.  Exiting.")
        return None
    #
    # Main column generation loop.   Optimize the restricted master,
    # query the dual values, create the subproblem, solve the subproblem
    # to obtain a new column to add to the restricted master.
    #
    quit    = False
    cgiters = 0
    cons   = rm.getConstrs()
    while not quit:
        if verbose:
            filename = "restrictedmaster" + str(cgiters) + ".lp"
            rm.write(filename)
        rm.optimize()
        if rm.status != GRB.OPTIMAL:
            print("Unexpected restricted master solve status ", rm.status)
            print("Exiting.")
            return None
        pivals = []
        for c in cons:
            pivals.append(c.PI)
        subprob = cg_subprob(stockwidth, cutwidths, pivals)
        if verbose:
            filename = "subprob" + str(cgiters) + ".lp"
            subprob.write(filename)

        subprob.optimize()
        if subprob.status != GRB.OPTIMAL:
            print("Unexpected restricted master solve status ", subprob.status)
            print("Colgen iteration ", cgiters, ". Exiting.")
            return None
        if subprob.objval <= 1.0:
            #
            # Restricted master is optimal.
            #
            print("**************************************************")
            print("Optimal solution found with objective ", rm.objval)
            print("Writing solution to file restrictedmaster.sol.")
            rm.write("restrictedmaster.sol")
            quit = True
        else:
            pattern = []
            vars    = subprob.getVars()
            for v in vars:
                pattern.append(v.X)
            #
            # We could let Gurobi's presolve derive the upperbound on the
            # variable associated with the new pattern, but we'll do it
            # here since the information is more directly available.
            bestub  = max(cutdemands)
            epint   = rm.getParamInfo("IntFeasTol")[2]
            newcol  = gp.Column()
            newcol.addTerms(pattern, cons)
            for i in range(numwidths):
                if pattern[i] <= epint:
                    continue
                newub = cutdemands[i]/pattern[i]
                if newub < bestub:
                    newub = bestub
            
            rm.addVar(lb=0.0, ub=bestub, obj=1.0, vtype=GRB.CONTINUOUS, \
                      name="pattern" + str(cgiters + 1 + numwidths), \
                      column=newcol)
            rm.update()
            if verbose:
                print("Pattern values added in CG Iteration ", cgiters,)
                for p in pattern:
                    print(p)
            cgiters += 1
    cgtime = time.time() - starttime
    print("Column Generation Time: ", cgtime)
    print("Column Generation Iterations: ", cgiters)    
    return rm
#
# Build the initial restricted master using the obvious feasible solution
# of assigning a single cutwidth to a single stock roll or board.  Can
# assign multiple rolls to a single stock if they fit (i.e. their width is
# at most half of the stock width).
#
def init_restricted_master(stockwidth, cutwidths, cutdemands):
    #
    # First determine the number of variables in the restricted master.
    #
    # import pdb; pdb.set_trace()
    rm        = gp.Model()
    numwidths = len(cutwidths)
    rm.addConstrs((gp.LinExpr() >= cutdemand for cutdemand in cutdemands), \
                  name="cutwidth")
    rm.update()
    cons      = rm.getConstrs()
    for i in range(numwidths):
        cutcol  = numwidths*[0]
        pattern = gp.Column()
        cwcount = math.floor(stockwidth/cutwidths[i])
        if cwcount == 0:
            print("Cut width ", i, "has length ", cutwidths[i], \
                  " which exceeds stock width of ", stockwidth)
            print("Original model is infeasible.  Exiting.")
            return None
        cutcol[i]  = cwcount
        pattern.addTerms(cutcol, cons)
        rm.addVar(lb=0.0, ub=math.ceil(cutdemands[i]/cwcount), obj=1.0, \
                  vtype=GRB.CONTINUOUS, name="pattern" + str(i+1), column=pattern)
    return(rm)

#
#   Column generation subproblem creation.
#   Given the dual variables pi, the stock width W, and the
#   cut widths wi, the subproblem is
#   max pi'z
#   s.t.
#   sum  zi*wi <= W
#    i
#   zi >= 0, integer
#
def cg_subprob(stockwidth, cutwidths, pivals):
    numvars = len(cutwidths)
    cgm     = gp.Model()
    cgmvars = cgm.addVars(numvars, lb=0, vtype=GRB.INTEGER, name="z")
    for i in range(numvars):
        cgmvars[i].obj = pivals[i]
        cgmvars[i].ub   = math.floor(stockwidth/cutwidths[i])
    lhs = gp.LinExpr()
    lhs.addTerms(cutwidths, cgmvars.values())
    cgm.addConstr(lhs <= stockwidth, name="knapsack")
    cgm.ModelSense = GRB.MAXIMIZE
    cgm.update()
    return cgm



#
#   Compact formulation:
#   xij = number of rolls of width cutwidths[i] cut from roll j
#   i   in {1,...,K}       where K = # of different cut widths
#   j   in {1,...,n}       where n is determined by cut widths and demands
#   yj  = 1 if roll j is used.
#         n
#   min  sum yj 
#        j=1
#   s.t.
#    K                     yj = 1 if any smaller rolls cut from stock roll j
#   sum  xij <= Uyj        j = 1,...,n     U determined by width data
#   i=1
#                          Sum of widths of smaller roll cut from stock roll j  
#    K                     cannot exceed the stock width
#   sum  wi*xij <= W       j = 1,...,n     w = cut widths, W = stock width
#   i=1 
#
#    n                     Total rolls of width i cut from stock >= demand
#   sum xij >= di          i = 1,...,K
#   j=1
#
#   xij >= 0, integer; yj binary
def compact_formulation(stockwidth, cutwidths, cutdemands, verbose=False):
    # import pdb; pdb.set_trace()
    #
    # First compute problem dimensions from data provided.
    # U, the max number of rolls that can be cut from a stock roll,
    # based cutting a single stock roll into rolls of the smallest
    # possible width
    #
    K = len(cutwidths)
    U = math.floor(stockwidth/min(cutwidths))
    #
    # n, the number of stock rolls from which cutting can be done,
    # is computed by considering the obvious feasible solution of
    # cutting only rolls of the same width from any stock roll.
    # 
    n = 0
    for i in range(K):
        num1 = math.floor(stockwidth/cutwidths[i])
        n   += math.floor(cutdemands[i]/num1)
    n = int(n)
    #
    # Ready to build model.  Create variables first.
    #
    m = gp.Model()
    y = m.addVars(n, lb=0, ub=1, vtype=GRB.BINARY, name="y")
    x = m.addVars(K, n, lb=0, ub=U, vtype=GRB.INTEGER, name="x")
    #
    # Create constraints and objective.
    #
    m.addConstrs((x.sum('*', j) <= U*y[j] for j in range(n)), name="rollused")
    for j in range(n):
        lhs = gp.LinExpr()
        lhs.addTerms(cutwidths, x.select('*',j))
        m.addConstr(lhs <= stockwidth, name = "widthlimit" + str(j))
        
    m.addConstrs((x.sum(i,'*') >= cutdemands[i] for i in range(K)), \
                 name="demand")
    m.setObjective(gp.quicksum(y))
                           
    
    return m
