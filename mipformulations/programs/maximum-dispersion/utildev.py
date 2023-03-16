#!/usr/bin/env python3.7

# Copyright 2020, Gurobi Optimization, LLC

#
# Display row or column nonzero count histograms for the specified model
# Type should be a string containing either "r", "c" or "rc" 
# to tell the program to provide histograms for rows, columns or both.
#

import gurobipy as gp
from gurobipy import GRB
from collections import Counter
from collections import defaultdict
from fractions import Fraction
import math
import random

CLIQUE  = 0
SETPACK = 1
SETPART = 2
SETCOV  = 3


def histogram(model, type="rc", plot=True):
    if not isinstance(type, str):
        print('Invalid histogram type input.  Should be "r", "c", "rc" or "cr".')
        print('Exiting')
        return False
              
    rcounts = []
    ccounts = []
#
#   Build list of nonzero counts for constraints and matrices.
#   Counter function creates a dictionary of the number of occurrences
#   of each nonzero count.
#  
    for i in range(len(type)):
        if type[i] == 'r':
            cons    = model.getConstrs()
            for c in cons:
                lhs = model.getRow(c)
                rcounts.append(lhs.size())
            rcounts.sort()
            rhist = Counter(rcounts)   # sparse dict of {value, count}
        elif type[i] == 'c':
            vars = model.getVars()
            for v in vars:
                col = model.getCol(v)
                ccounts.append(col.size())
            ccounts.sort()
            chist = Counter(ccounts)   # sparse dict of {value, count}
#
#   Display the histogram to the console.  Includes some feeble attempts
#   to get consistent alignment in the output based on the numeric values
#   to print.  Do the rows first.
#
    if len(rcounts) > 0:
        print ("Histogram of row nonzero counts:")
        maxrlen    = max(rcounts)
        maxrcount  = max(rhist.values())
        width      = 2 + math.ceil(math.log(float(max([maxrlen, maxrcount])),
                                            10.0))
        numperline = int (64/width)
        itemcount  = 0
        line1 = "Nonzero Count:  "
        line2 = "Number of Rows: "
        for nzcnt in rhist.keys():
            line1 += (str(nzcnt)).rjust(width, ' ')
            line2 += (str(rhist[nzcnt])).rjust(width, ' ')
            itemcount += 1
            if itemcount == numperline:
                print(line1)
                print(line2)
                itemcount = 0
                line1 = "Nonzero Count:  "
                line2 = "Number of Rows: "
        if itemcount > 0:     # last line may have fewer elements to print
            print(line1)
            print(line2)
            
#
#   Rows done; do likewise with columns.  TODO:  code is similar enough
#   that it should be integrated into a single function call that will
#   do either rows or columns.   Probably just need to pass in 3 or 4
#   arguments to the function.
#

    if len(ccounts) > 0:
        print ("Histogram of column nonzero counts:")
        maxclen    = max(ccounts)
        maxccount  = max(chist.values())
        width      = 2 + math.ceil(math.log(float(max([maxclen, maxccount])),
                                            10.0))
        numperline = int (64/width)
        itemcount  = 0
        line1 = "Nonzero Count:  "
        line2 = "Number of Cols: "
        for nzcnt in chist.keys():
            line1 += (str(nzcnt)).rjust(width, ' ')
            line2 += (str(chist[nzcnt])).rjust(width, ' ')
            itemcount += 1
            if itemcount == numperline:
                print(line1)
                print(line2)
                itemcount = 0
                line1 = "Nonzero Count:  "
                line2 = "Number of Cols: "
        if itemcount > 0:     # last line may have fewer elements to print
            print(line1)
            print(line2)

#
#   Now plot the histograms.  One bin for each distinct nonzero count
#

    if plot == False:
        return
    import pdb; pdb.set_trace()

    try:
        import pandas as pd
    except ImportError as e:
        print("Need pandas installed to plot histogram.")
        return
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:
        print("Need matplotlib.pyplot installed to plot histogram.")
        return
        

    
    if len(type) == 2:
        fig, axes = plt.subplots(nrows=2)
    else:
        fig, axes = plt.subplots(nrows=1)
    if len(rcounts) > 0:
        df      = pd.DataFrame(rcounts, columns=['Row Nonzero Counts'])
        if len(type) == 1:
            rowhist = df.plot.hist(bins=len(rhist), ax=axes)
        else:
            rowhist = df.plot.hist(bins=len(rhist), ax=axes[0])

    if len(ccounts) > 0:
        df      = pd.DataFrame(ccounts, columns=['Col Nonzero Counts'])
        if len(type) == 1:
            colhist = df.plot.hist(bins=len(chist), ax=axes)
        else:
            colhist = df.plot.hist(bins=len(chist), ax=axes[1])
            
    plt.show()



#
#   Translate a single objective model into a multi objective model
#   based on the list in the thresholds argument.  Letting t1,...,tk
#   be the k threshold values, we create one objective consisting of
#   all variables with obj coeffs <= t1.   The next objective consists of
#   all variables with obj coeffs in (t1, t2] and so on.   The objective
#   with largest coefficients gets solved first.   The rescale parameter
#   indicates whether to rescale the coefficients.   This is a benefit
#   of multi obj solves, as it can help avoid numerical troubles due to
#   large objective coefficients (either by themselves or mised with small
#   ones.   The default of 1 scales by the maximum element, 2 scales by
#   the geometric mean of the coefficients associated each multi obj
#   solve, and setting rescale to 0 leaves the objective unscaled.
#
    
#
#   Check with Kostja; Tobi says this should already 
#def grbtomulti(model, thresholds, rescale=1):


#
#   Takes a linear combination of the constraints listed in cons
#   mults is an optional list of same length as cons containing
#   the multipliers that specify the linear combination.  If no
#   multipliers specified, they are all assumed to be 1.0
#   Returns a tuple consisting of aggregated left hand side, aggregated
#   sense, aggregated right hand side
#
#   If a file is specified, is assumed to contain a list of constraint
#   names and optional multipliers, in which case we use that input
#   to combine the constraints.
#
def combineconstraints(model, mcons, mmults=None, filepath = None, showit=True):
    
    if (filepath):
        with open(filepath) as file:
            condict = {}
            for c in model.getConstrs():
                condict[c.ConstrName] = c

            fcons  = []
            fmults = []
            for i, line in enumerate(file):
                tokens = line.split()
                if len(tokens) < 1 or len(tokens) > 2:
                    print("Warning: bogus line ", i, ", ignoring")
                    continue
                fcons.append(condict[tokens[0]])
                if len(tokens) == 2:
                    fmults.append(float(tokens[1]))
                else:
                    fmults.append(1.0)
            cons  = fcons
            mults = fmults
    else:
        cons  = mcons
        mults = mmults
                
    baselhs    = gp.LinExpr()
    baserhs    = 0
    firstsense = cons[0].sense
    showsense  = True
    if mults == None:            # sum up the constraints
        for c in cons:
            
            baselhs.add(model.getRow(c))
            baserhs += c.rhs
            if c.sense != firstsense:
                showsense = False
    else:                        # more general linear combination
        for c,t in zip(cons, mults):
            baselhs.add(model.getRow(c), t)
            baserhs += t*c.rhs
            if c.sense != firstsense:
                showsense = False
#
#   Per Gurobi docs, left hand side doesn't yet combine terms containing
#   the same variables; that happens when left hand side is added as
#   a constraint.   For example, if we added constraints x1 + x2 <= 1 and
#   x2 + x3 <= 1, the baselhs at this point would be x1 + x2 + x2 + x3
#   rather than x1 + 2x2 + x3.  Make the consolidation of terms happen now
#   in a separate, dummy constraint
#

    model.addConstr(baselhs, GRB.EQUAL, baserhs)
    model.update()
    cons = model.getConstrs()
    baselhs = model.getRow(cons[len(cons) - 1])  # Consolidation of terms
#
#   For now we only display the constraint sense if all senses in
#   constraints to be combined are the same.   TODO: Could extend to
#   provide sense as long as we don't mix <= and >= constraints.
#
                
    if showsense:
        sense = str(firstsense) + "="
    else:
        sense = "~~"      # mix of constraints
        
    if showit:
        print ("Combined Constraint: ")
        print (baselhs, "  ", sense, "  ", str(baserhs))

    model.remove(cons[len(cons) - 1])       # Restore model
    model.update()
    return (baselhs, sense, baserhs)


#
#   Extract linear constraints with names specified by prefix and optional
#   suffix
#

def extractconstraints(model, prefix, suffix=None):

    cons = []
    for c in model.getConstrs():
        cstr = c.ConstrName
        if prefix != None:
            if not cstr.startswith(prefix):
                continue
        if suffix != None:
            if not cstr.endswith(suffix):
                continue
        cons.append(c)

    return cons

#
#   Extract quadratic constraints with names specified by prefix
#   and optional suffix
#

def extractqconstraints(model, prefix, suffix=None):

    qcons = []
    for qc in model.getQConstrs():
        cstr = qc.QCName
        if prefix != None:
            if not cstr.startswith(prefix):
                continue
        if suffix != None:
            if not cstr.endswith(suffix):
                continue
        qcons.append(qc)

    return qcons


#
#   Extract variables with names specified by prefix and optional suffix
#   Also allow for extracting by upper or lower bounds above/below a
#   specified threshold.  Also allow for extracting by a string of types
#   consisting of a subset of "BIC" (Binary/Integer/Continuous)
#

def extractvariables(model, prefix=None, suffix=None, bndthresh=None, \
                     types=None):

    vars = []
    for v in model.getVars():
        vstr = v.VarName
        if prefix != None:
            if not vstr.startswith(prefix):
                continue
        if suffix != None:
            if not vstr.endswith(suffix):
                continue
        if bndthresh != None:
            if (v.UB >= GRB.INFINITY or v.UB < bndthresh) and  \
               (v.LB <= -GRB.INFINITY or v.LB > -bndthresh) :
                continue
        if types != None:
            hit = False
            for char in types:
                if v.Vtype == char:
                    hit = True
                    break
            if not hit:
                continue
        vars.append(v)

    return vars
#
#   Returns the variables contained in a linear expression
#
def extractvariablesfromexpr(linexpr):
    retvars = []
    for j in range(linexpr.size()):
        retvars.append(linexpr.getVar(j))

    return retvars
    

#
#   Display the constraints intersecting the specified variable.
#

def displayvariable(model, var):
    
    col = model.getCol(var)
    print("Constraints intersecting variable ", var.VarName, ":")
    for i in range(col.size()):
        displayconstraint(model, col.getConstr(i))


def displayconstraint(model, con):
    print(con.ConstrName, ": ", model.getRow(con), con.sense, con.rhs)

    
        
def buildvardict(model):
    vardict = {}
    for v in model.getVars():
        vardict[v.VarName] = v
    return vardict
#
#   Gives priority specified in pvalue to all variables in model
#   that start with the specified prefix and end with the optional
#   specified suffix.  Either prefix or suffix can be None   
#
def setpriorities(model, prefix, pvalue, suffix = None, write = False):
    for v in model.getVars():
        pstr = v.VarName
        if prefix != None:
            if not pstr.startswith(prefix):
                continue
        if suffix != None:
            if not pstr.endswith(suffix):
                continue
        if v.Vtype == GRB.CONTINUOUS:       # No branching on cont. vars
            print("Priority assigned to continuous variable " + v.VarName)
            print("Ignoring.")
        v.BranchPriority = pvalue

    model.update()
        
    if write:
        if len(model.ModelName) == 0:
            modelname = "model"
        else:
            modelname = model.ModelName
        model.write(modelname + "_pri.attr")


#
#   Alternate priority setting scheme.   Extract a string based on
#   the start and end position (inclusive) and use that as the
#   priority for the variable.  
#
def setstringpriorities(model, startpos, endpos, reverse=False, write=False):
    import pdb; pdb.set_trace()
    maxpri = 0
    for v in model.getVars():
        if v.Vtype == GRB.CONTINUOUS:       # No branching on cont. vars
            continue
        if endpos > len(v.VarName):
            continue
        pstr   = v.VarName[startpos: endpos + 1]
        pvalue = int(pstr)
        if pvalue > maxpri:
            maxpri = pvalue
        v.BranchPriority = int(pvalue)

    #
    # Allow for giving highest priority to lowest value in the extracted
    # strings (e.g., you want to give highest priority to time periods
    # with the lowest values.
    #
    if reverse:
        model.update()
        for v in model.getVars(): 
            if v.Vtype == GRB.CONTINUOUS:    
                continue
            v.BranchPriority = maxpri - v.BranchPriority
           
    model.update()
        
    if write:
        if len(model.ModelName) == 0:
            modelname = "model"
        else:
            modelname = model.ModelName
        model.write(modelname + "_pri.attr")

        
#
#   Extract the binary variables from the model.
#
def getbinvars(model):
    binvars = []
    for v in model.getVars():
        type = v.Vtype
        if v.Vtype == GRB.CONTINUOUS:
            continue
        elif v.Vtype == GRB.BINARY or (v.Vtype == GRB.INTEGER and v.UB == 1):
            binvars.append(v)
    return binvars
        
#
#   Gives higher priorities to binaries.  pvalue specifies
#   priorities for binaries (all other discrete objects have
#   priority 0).  
#
def binvarsfirst(model, pvalue, write = False):
    for v in model.getVars():
        type = v.Vtype
        if v.Vtype == GRB.CONTINUOUS:
            continue
        elif v.Vtype == GRB.BINARY or (v.Vtype == GRB.INTEGER and v.UB == 1):
            v.BranchPriority = pvalue

    model.update()
        
    if write:
        if len(model.ModelName) == 0:
            modelname = "model"
        else:
            modelname = model.ModelName
        model.write(modelname + "_pri.attr")
        

#
#   Extract the general integer variables from the model.
#
def getintvars(model):
    intvars = []
    for v in model.getVars():
        type = v.Vtype
        if v.Vtype == GRB.CONTINUOUS:
            continue
        elif v.Vtype == GRB.INTEGER:
            intvars.append(v)
    return intvars


#
#   Gives higher priorities to variables with nonzero objective coefficients.
#   pvalue specifies priorities for variables with objective coeffs.
#   All other variables get priority 0.  
#
def objvarsfirst(model, pvalue, write = False):
    for v in model.getVars():
        type = v.Vtype
        if v.Vtype == GRB.CONTINUOUS:
            continue
        elif v.Vtype == GRB.BINARY or (v.Vtype == GRB.INTEGER and v.Obj != 0.0):
            v.BranchPriority = pvalue

    model.update()
        
    if write:
        if len(model.ModelName) == 0:
            modelname = "model"
        else:
            modelname = model.ModelName
        model.write(modelname + "_pri.attr")


#
#   Scales the objective in the model m by the factor specified
#   in mult.   Alters the model; use on a copy if you want to
#   preserve the original model.  
#
def scaleobj(model, mult):
    objexpr = model.getObjective()
    linpart = None
    if isinstance(objexpr, gp.LinExpr):
        linpart = objexpr
    elif isinstance(objexpr, gp.QuadExpr):
        linpart = objexpr.getLinExpr()
    else:
        print("Objective is neither quadratic nor linear.  No scaling done.")
        return None
    scaledlinpart  = None
    if linpart.size() > 0:
        scaledlinpart  = scalelinexpr(linpart, mult)
    scaledquadpart = None
    if isinstance(objexpr, gp.QuadExpr):
        scaledquadpart = scalequadexpr(objexpr, mult)
#
#   Scaling done; update the objective in the model
#
    if scaledquadpart == None:
        model.setObjective(scaledlinpart)
    else:
        scaledquadexpr = gp.QuadExpr()
        scaledquadexpr.add(scaledquadpart)
        if scaledlinpart != None:
            scaledquadexpr.add(scaledlinpart)
        model.setObjective(scaledquadexpr)
    model.update()
    scaledobjexpr = model.getObjective()
    return scaledobjexpr

def scaleconstraint(model, con, mult):
    lhs         = model.getRow(con)
    rhs         = con.rhs
    scaledlhs   = scalelinexpr(lhs, mult)
    scaledrhs   = rhs*mult
    scaledsense = con.sense
    scaledcon   = None
    if mult >= 0:
        scaledcon = (scaledlhs, scaledsense, scaledrhs)
    else:
        if con.sense == '<':
            scaledsense = '>'
        elif con.sense == '>':
            scaledsense = '<'            
        scaledcon = (scaledlhs, scaledsense, scaledrhs)

    return scaledcon
        
    

def scalequadexpr(quadexpr, mult):
    scaledqexpr = gp.QuadExpr()
    for k in range(quadexpr.size()):
        qcoef = mult*quadexpr.getCoeff(k)
        xi    = quadexpr.getVar1(k)
        xj    = quadexpr.getVar2(k)
        scaledqexpr += qcoef*xi*xj
    return scaledqexpr

def scalelinexpr(linexpr, mult):
    scaledexpr = gp.LinExpr()
    for j in range(linexpr.size()):
        var   = linexpr.getVar(j)
        coeff = linexpr.getCoeff(j)
        coeff *= mult
        scaledexpr.addTerms(coeff, var)
    return scaledexpr

#
#   Zeros out elements of a vector with absolute value.  Alters the
#   input vector; use a copy if you want to preserve the original
#   vector.
#
def zerovec(vec, zapval):
    for i in range(len(vec)):
        if vec[i] < zapval:
            vec[i] = 0

#
#   Zeros out objective coefficients in the model m by the threshold value
#   specified in zapval.   Alters the model; use on a copy if you want to
#   preserve the original model.  
#
def zeroobj(model, zapval):
    objexpr = model.getObjective()
    linpart = None
    if isinstance(objexpr, gp.LinExpr):
        linpart = objexpr
    elif isinstance(objexpr, gp.QuadExpr):
        linpart = objexpr.getLinExpr()
    else:
        print("Objective is neither quadratic nor linear.  No scaling done.")
        return None
    zeroedlinpart  = None
    if linpart.size() > 0:
        zeroedlinpart  = zerolinexpr(linpart, zapval)
    zeroedquadpart = None
    if isinstance(objexpr, gp.QuadExpr):
        zeroedquadpart = zeroquadexpr(objexpr, zapval)
#
#   Scaling done; update the objective in the model
#
    if zeroedquadpart == None:
        model.setObjective(zeroedlinpart)
    else:
        zeroedquadexpr = gp.QuadExpr()
        zeroedquadexpr.add(zeroedquadpart)
        if not (zeroedlinpart == None):
            zeroedquadexpr.add(zeroedlinpart)
        model.setObjective(zeroedquadexpr)
    model.update()
    zeroedobjexpr = model.getObjective()
    return zeroedobjexpr

def zeroquadexpr(quadexpr, zapval):
    zeroedqexpr = gp.QuadExpr()
    for k in range(quadexpr.size()):
        qcoef = quadexpr.getCoeff(k)
        if abs(qcoef) >= zapval:
            xi    = quadexpr.getVar1(k)
            xj    = quadexpr.getVar2(k)
            zeroedqexpr += qcoef*xi*xj
    return zeroedqexpr

def zerolinexpr(linexpr, zapval):
    zeroedexpr = gp.LinExpr()
    for j in range(linexpr.size()):
        coeff = linexpr.getCoeff(j)
        if abs(coeff) >= zapval:
            var   = linexpr.getVar(j)
            zeroedexpr.addTerms(coeff, var)
    return zeroedexpr




#
#   Perturbs the objective in the model by the amount specified
#   in pert multiplied by a random value in(0,1).   Alters the model;
#   use on a copy if you want to preserve the original model.  
#
def perturbobj(model, pert=1e-04):
    objexpr = model.getObjective()
    linpart = None
    if isinstance(objexpr, gp.LinExpr):
        linpart = objexpr
    elif isinstance(objexpr, gp.QuadExpr):
        linpart = objexpr.getLinExpr()
    else:
        print("Objective is neither quadratic nor linear.  No scaling done.")
        return None
    perturbedlinpart  = None
    if linpart.size() > 0:
        perturbedlinpart  = perturblinexpr(linpart, pert)
    perturbedquadpart = None
    if isinstance(objexpr, gp.QuadExpr):
        perturbedquadpart = perturbquadexpr(objexpr, pert)
#
#   Perturbation done; update the objective in the model
#
    if perturbedquadpart == None:
        model.setObjective(perturbedlinpart)
    else:
        perturbedquadexpr = gp.QuadExpr()
        perturbedquadexpr.add(perturbedquadpart)
        if perturbedlinpart != None:
            perturbedquadexpr.add(perturbedlinpart)
        model.setObjective(perturbedquadexpr)
    model.update()
    perturbedobjexpr = model.getObjective()
    return perturbedobjexpr

def perturbquadexpr(quadexpr, pert):
    perturbedqexpr = gp.QuadExpr()
    for k in range(quadexpr.size()):
        qcoef = quadexpr.getCoeff(k) + pert*random.uniform(0.0, 1.0)
        xi    = quadexpr.getVar1(k)
        xj    = quadexpr.getVar2(k)
        perturbedqexpr += qcoef*xi*xj
    return perturbedqexpr

def perturblinexpr(linexpr, pert):
    perturbedexpr = gp.LinExpr()
    for j in range(linexpr.size()):
        var   = linexpr.getVar(j)
        coeff = linexpr.getCoeff(j)
        coeff += pert*random.uniform(0.0, 1.0)
        perturbedexpr.addTerms(coeff, var)
    return perturbedexpr




#
#   Checks if constraint is a set pack/partition/cover constraint, i.e.
#   sums of binaries <= 1 or == 1 or >= 1.  By default, looks for CLIQUE
#   constraints, i.e. either packing or partitioning.  Accepts hidden
#   versions of set p/p/c constraints where all coefficients are -1 and
#   the rhs is -1.
#
def issetconstraint(model, con, type=CLIQUE):
    sense    = con.Sense
    rhs      = con.RHS
    lhs      = model.getRow(con)
    coeffs   = []
    maxcoeff = -math.inf
    mincoeff = math.inf
    allbin   = True
    for j in range(lhs.size()):
        coeff = lhs.getCoeff(j)
        if coeff > maxcoeff:
            maxcoeff = coeff
        if coeff < mincoeff:
            mincoeff = coeff
        coeffs.append(coeff)
        var   = lhs.getVar(j)
        isbin = var.Vtype == GRB.BINARY
        if var.Vtype == GRB.INTEGER:
            isbin = False
            if var.LB == 0 and var.UB == 1:
                isbin = True
        allbin = allbin and isbin    
                
    if abs(maxcoeff - mincoeff) > 0.0:  
        return False
    if maxcoeff != rhs:
        return False
    if not allbin:
        return False
    if type == CLIQUE or type == SETPART:
        if maxcoeff < 0 and sense == '<':
            return False
        if maxcoeff > 0 and sense == '>':
            return False
    elif type == SETPACK:
        if maxcoeff < 0 and sense == '<':
            return False
    elif type == SETCOV:
        if maxcoeff < 0 and sense == '>':
            return False        
#
#   All our filters have been passed.   we have a constraint with all
#   all binary variables, a single nonzero coefficient value for all lhs
#   terms and the right hand side.   And the constraint sense ensures
#   that we have a clique constraint
#
    return True





#
#   Translates any clique constraints to SOS constraints.  Alters model;
#   create a copy if you want to preserve the original model.  By default
#   will replace the clique constraints with SOS constraints; set delete
#   to False if you want to retain the clique the clique constraints in
#   the model.
#
def cliquestoSOS(model, delete=True):
    cliquecons = []
    for c in model.getConstrs():           # extract the clique constraints
        if iscliqueconstraint(model, c):
            cliquecons.append(c)
#
#   For each clique constraint, extract the variables into the list
#   that identifies the variables in the SOS.   Use standard weights
#   1, 2, 3,...
#
    for cl in cliquecons:
        lhs = model.getRow(cl)
        varlist    = []
        weightlist = []
        for j in range(lhs.size()):
            varlist.append(lhs.getVar(j))
            weightlist.append(j+1)
        model.addSOS(GRB.SOS_TYPE1, varlist, weightlist)
    if delete:
        model.remove(cliquecons)
        
            





#
#   Looks for all binary variable knapsack constraints.  Does not
#   consider complements; only checks if original variables form
#   a knapsack constraint.
#
def isknapsackconstraint(model, con):
    sense    = con.Sense
    rhs      = con.RHS
    lhs      = m.getRow(con)
    coeffs   = []
    maxcoeff = -math.inf
    mincoeff = math.inf
    allbin   = True
    for j in range(lhs.size()):
        coeff = lhs.getCoeff(j)
        if coeff > maxcoeff:
            maxcoeff = coeff
        if coeff < mincoeff:
            mincoeff = coeff
        coeffs.append(coeff)
        var   = lhs.getVar(j)
        isbin = var.Vtype == GRB.BINARY
        if var.Vtype == GRB.INTEGER:
            isbin = False
            if var.LB == 0 and var.UB == 1:
                isbin = True
        allbin = allbin and isbin    
                
    if not allbin:
        return False
    if maxcoeff < 0 and sense == '<':
        return False
    if maxcoeff > 0 and sense == '>':
        return False
#
#   All our filters have been passed.   we have a constraint with all
#   all binary variables, a single nonzero coefficient value for all lhs
#   terms and the right hand side.   And the constraint sense ensures
#   that we have a clique constraint
#
    return True



        
#
#   More detailed statistics beyond the printStats() function
#
def moreStats(model):
#
#   Linear constraints by type
#    
    cons = model.getConstrs()
    contypes = [0, 0, 0]    # respective number of < , > and = constraints
    typedict = {'<' : 0, '>' : 1, '=' : 2}
    for c in cons:
        contypes[typedict[c.Sense]] += 1
    nle  = contypes[0]
    nge  = contypes[1]
    neq  = contypes[2]
    ntot = sum(contypes)
#
#   Variables by type.
#
        
#
#   Results
#
    print("Additional statistics for model ", model.ModelName)
    print("Linear Constraints: ", nle, " <=; ", neq, " ==; ", nge, " >=.")
    print(ntot, " total linear constraints.")


#
#   Print all nonzero solution values, optionally specifying an absolute
#   threshold only above which the value is printed.
#
def printx(m, thresh=0.0, varstoprint=None):
    if varstoprint == None:
        varstoprint = m.getVars()
        
    for v in varstoprint:
        if abs(v.X) > thresh:
            print("Variable ", v.varName, " has value ", v.X)
            
#
#   Print all nonzero dual values, optionally specifying an absolute
#   threshold only above which the value is printed.
#
def printpi(m, thresh=0.0, constoprint=None):
    if constoprint == None:
        cons = m.getConstrs()
    else:
        cons = constoprint
        
    if hasattr(cons[0], 'CBasis'): 
        for c in cons:
            if abs(c.PI) > thresh:
                print("Constraint ", c.constrName, " has value ", c.Pi)
    else:
        print("No Basis currently available.")
        
#
#   Print all nonzero solution values, optionally specifying an absolute
#   threshold only above which the value is printed.
#
def printrc(m, thresh=0.0):
    vars = m.getVars()
    if hasattr(vars[0], 'VBasis'): 
        for v in vars:
            if abs(v.RC) > thresh:
                print("Variable ", v.varName, " has reduced cost ", v.RC)
    else:
        print("No Basis currently available.")

#
#   Count the number of nonzeros in the objective (nonzero tolerance
#   optional
#
def countobjnonzeros(model, threshold=0.0):
    vars = model.getVars()
    count = 0
    for v in vars:
        if abs(v.Obj) > threshold:
            count += 1
    return count
        
#
#   Checks for possible penalty variables (column singletons with
#   nonzero objective coefficients.  Does not look at coefficient sign
#   or constraint sense to verify that the candidate is truly a penalty
#   variable.
#
def checkforpenaltyvars(model):
    vars = model.getVars()
    for v in vars:
        col = model.getCol(v)
        if col.size() == 1 and v.Obj == 0.0:
            print(v.VarName, " appears to be a penalty variable.")
#
#   Simple dot product between two vectors
#       
def dot(l1, l2):
    if len(l1) != len(l2):
        print("Input arrays must have same length")
        return 0
           
    result = sum(l1[i]*l2[i] for i in range(len(l1)))
    return result

#
#   Compute distance between points.  Default norm is L2, but supports
#   L0, L1 and LInf as well.  p1 and p2 are lists of values defining the
#   points.
#
def distance(p1, p2, norm="L2", tol=None):
    n1 = len(p1)
    if n1 != len(p2):
        print("Mismatch; points must have the same length. ")
        print("p1 had length ", n1, "; p2 has length ", len(p2))
        return 0.0
    dist = 0.0
    if norm == "L2":
        for j in range(n1):
            dist += (p1[j] - p2[j])**2
        dist = math.sqrt(dist)
    elif norm == "L1":
        for j in range(n1):
            dist += abs(p1[j] - p2[j])
    elif norm == "L0":
        if tol==None:
#
#       Calculate the zero tolerance based on the magnitudes of the
#       input data.  Assumes 64 bit doubles and hence 16 digits of
#       accuracy.
#
            avg    = (sum(p1) + sum(p2)) / n1
            digits = math.ceil(math.log(avg, 10))
            tol    = math.pow(10, (digits - 16))
        count = 0
        for j in range(n1):
            if abs(p1[j] - p2[j]) < tol:
                count += 1
            dist = count
    elif norm == "LInf":
        max = 0.0
        for j in range(n1):
            t = abs(p1[j] - p2[j])
            if t > max:
                max = t
        dist = max
    else:
        print("Unknown norm specified. Acceptable norms are L0, L1, L2, LInf")
        print("You specified ", norm)
        return 0
    
    return dist

#
#  Convert CPLEX XML based MST file to Gurobi .sol file format
#        
def cpxmst2grbsol(mstinfilename, soloutfilename):
    try:
        soloutfile = open(soloutfilename, "w")
    except IOError:
        print("Unable to open ", soloutfilename, " for writing")
        return
    
    try:
        with open(mstinfilename) as file:
            for i, line in enumerate(file):
                tokens = line.split()
                if tokens[0].startswith("problemName"): # e.g. problemName="zzz.mps"
                    kstart   = tokens[0].find("=\"") + 2
                    kend     = tokens[0][kstart:].find("\"")
                    probname = tokens[0][kstart:kstart + kend]
                    soloutfile.write("# Solution for model " + probname + "\n")
                elif tokens[0] == "<variable":
                    kstart   = tokens[1].find("\"") + 1
                    kend     = tokens[1][kstart:].find("\"")
                    varname  = tokens[1][kstart:kstart + kend]
                    kstart   = tokens[3].find("\"") + 1
                    kend     = tokens[3][kstart:].find("\"")
                    varval   = tokens[3][kstart:kstart + kend]
                    soloutfile.write(varname +  " " + varval + "\n")
                elif tokens[0] == "</variables>":   # Done with variables
                    soloutfile.close
                    break                          # Ignore rest of file
    except IOError:
        print("File ", mstinfilename, " not found") 

#
#  This adds a constraint on the objective function to the model
#  This will be a <= constraint for minimization problems and a >=
#  constraint for maximization problems.  The most common use case
#  is to create an infeasible model and then generate the IIS.
#  The compute argument will do this by default
#
def conobj(model, conobjval, compute=True):
    objexpr  = model.getObjective()
    consense = '<'
    if model.ModelSense == GRB.MAXIMIZE:
        consense = '>'
    model.addConstr(objexpr, consense, conobjval, "_conobj")
    filename = model.ModelName + "_conobj"
    model.write(filename + ".lp")
    if compute:
        model.computeIIS()
        model.write(filename + ".ilp")
#
#   From Kostja's grbmodreport.py file in his github.
#        
def get_obj_c_frequencies(model, basis):

    int_count = defaultdict(int)
    cont_count = defaultdict(int)

    for x in model.getVars():
        if x.Obj != 0:
            if x.VType == 'C':
                cont_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1
            else:
                int_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1

    if len(cont_count) > 0:
        cont_result = [ [exponent, cont_count[exponent]] for exponent in range(min(cont_count.keys()), max(cont_count.keys()) + 1) ]
    else:
        cont_result = None
    if len(int_count) > 0:
        int_result = [ [exponent, int_count[exponent]] for exponent in range(min(int_count.keys()), max(int_count.keys()) + 1) ]
    else:
        int_result = None

    return cont_result, int_result
    


if __name__ == "__main__":
    # Use as customized command line tool
    import sys
    m  = gp.read(sys.argv[1])

#
#   Checks to see if the cuts in cutlist are globally valid.   Reverses
#   the constraint sense, tightens the right hand side by rhsdelta (either
#   provided by user or computed based on the rhs value of the candidate
#   cut, then confirms infeasibility.  For inequality cuts only.
#   NOTE: this will modify the model, so should run this on a copy of
#   the original model if you wish to preserve it.
#
def checkcuts(m, cutlist, rhsdelta=None, TimeLim=300.0):
    m.setParam("SolutionLimit", 1)    
    m.setParam("TimeLimit", TimeLim)
    numfails = 0
    for con in cutlist:
        oldrhs   = 0
        oldsense = None
        conname  = None
        isQC     = isinstance(con, gp.QConstr)
        if isQC:
            oldrhs   = con.QCRHS
            oldsense = con.qcsense
            conname  = con.QCName
        else:
            oldrhs   = con.RHS
            oldsense = con.sense
            conname  = con.ConstrName
            
        conrhs   = oldrhs
        consense = oldsense
             
        if rhsdelta == None:   # Need to compute RHS for reversed inequality
            if conrhs.is_integer():
                rhsdelta = 1.0
            else:
                rhsdelta = abs(conrhs)*0.1
                
        if consense == GRB.LESS_EQUAL:
            if isQC:
                con.qcsense = GRB.GREATER_EQUAL
                con.QCRHS += rhsdelta
            else:
                con.sense = GRB.GREATER_EQUAL
                con.rhs   += rhsdelta
        elif consense == GRB.GREATER_EQUAL:
            if isQC:
                con.qcsense = GRB.LESS_EQUAL
                con.QCRHS  -= rhsdelta
            else:
                con.sense = GRB.LESS_EQUAL
                con.RHS  -= rhsdelta
                
        else:
            print("Equality constraint ", con.ConstName, ": ignoring.")

        m.optimize()
        if m.Status != GRB.INFEASIBLE and m.Status != GRB.INF_OR_UNBD and \
           m.Status != GRB.TIME_LIMIT:
            print("Unexpected Optim. Status for ", conname, ": ", \
                  m.Status)
            print("Continuing.")
            numfails += 1
#            if numfails == 1:
#                m.write("qtest.lp")
        else:
            print("Cut for ", conname, " valid.")

#
#       Restore reversed constraint.
#
        if isQC:
            con.QCRHS   = oldrhs
            con.qcsense = oldsense
        else:
            con.RHS   = oldrhs
            con.sense = oldsense

    print ("Total number of invalid cuts: ", numfails)
#
#   From Kostja's model analyzer to get histogram of model frequencies
#
def get_obj_c_frequencies(model, basis):
    int_count = defaultdict(int)
    cont_count = defaultdict(int)
    int_result = []
    cont_result = []

    for x in model.getVars():
        if x.Obj != 0:
            if x.VType == 'C':
                cont_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1
            else:
                int_count[int(math.floor(math.log(abs(x.Obj), basis)))] += 1

    if len(cont_count) > 0:
        cont_result = [ [exponent, cont_count[exponent]] for exponent in range(min(cont_count.keys()), max(cont_count.keys()) + 1) ]

    if len(int_count) > 0:
        int_result = [ [exponent, int_count[exponent]] for exponent in range(min(int_count.keys()), max(int_count.keys()) + 1) ]

    return cont_result, int_result

def extractAmat(model):
    Amat = []
    for con in model.getConstrs():
        rowlhs = model.getRow(con)
        for j in range(rowlhs.size()):
            Amat.append(rowlhs.getCoeff(j))
    return Amat

#
#   Generic version of Kostja's obj_frequencies routine for an
#   arbitrary vector

def get_vector_frequencies(vec, base=10):
    count = defaultdict(int)
    result = []

    for v in vec:
        if v != 0:
            count[int(math.floor(math.log(abs(v), base)))] += 1

    if len(count) > 0:
        result = [ [exponent, count[exponent]] for exponent in \
                   range(min(count.keys()), max(count.keys()) + 1) ]

    return result
#
#   Prints out the number of different general constraint types from
#   the supplied list of general constraints (e.g. pass model.getGenConstrs()
#   to get info on all general constraints in the model).
#   Also returns a dictionary of the counts for each general constraint type.
#
def genConstrinfo(gcs):
    gctypedict = {GRB.GENCONSTR_MAX: "Max:", GRB.GENCONSTR_MIN: "Min:", \
                  GRB.GENCONSTR_ABS: "Abs:", GRB.GENCONSTR_AND: "And:", \
                  GRB.GENCONSTR_OR: "Or:", \
                  GRB.GENCONSTR_INDICATOR: "Indicator:",\
                  GRB.GENCONSTR_PWL: "Piecewise Linear:", \
                  GRB.GENCONSTR_POLY: "Polynomial:", \
                  GRB.GENCONSTR_EXP: "e^x:", \
                  GRB.GENCONSTR_EXPA: "a^x:",  GRB.GENCONSTR_LOG: "log_e(x):", \
                  GRB.GENCONSTR_LOGA: "log_a(x):", GRB.GENCONSTR_POW: "x^a:", \
                  GRB.GENCONSTR_SIN: "Sine:", GRB.GENCONSTR_COS: "Cosine:", \
                  GRB.GENCONSTR_TAN: "Tangent"}
    typecountdict = {}
    gctypes       = gctypedict.keys()
    maxstrlen     = 0
    for type in gctypes:
        if len(gctypedict[type]) > maxstrlen:
            maxstrlen = len(gctypedict[type])
        typecountdict[type] = 0
               
    for gc in gcs:
        typecountdict[gc.GenConstrType] += 1
        
    print("Total General Constraints: ", len(gcs))
    for type in gctypes:
        if typecountdict[type] > 0:
            print(gctypedict[type].ljust(maxstrlen, ' '), typecountdict[type])
            
    return typecountdict
#
#   Runs numruns different seeds, starting with seed startseed.  Leaves
#   original model untouched.
#
def runseeds(model, startseed, numruns, prmfile=None, attrfile=None,
             logfilename=None):
    m = model.copy()
    if prmfile != None:
        m.read(prmfile)
    for j in range(numruns):
        m.reset()
        if attrfile != None:
            m.read(attrfile)
        logname = "LogFile"
        if logfilename == None:
            logname += m.ModelName
        else:
            logname += logfilename
        m.setParam("LogFile", logname +  "_Seed" + str(startseed+j) + ".log")
        m.setParam("Seed", startseed + j)
        m.optimize()
#
#   Performs a single, constraint based, sifting solve.   In other
#   words, specify a subset of constraints to optimize first, then
#   add the remaining constraints and solve with dual simplex, since
#   we have a dual feasible starting basis.  Subset of constraints can
#   be specified as a list of constraints in mcons, or by a file constaining
#   constraint names.
#
def onesiftsolve(model, mcons, filepath=None):
#    import pdb; pdb.set_trace()
    #
    # Determine subset of constraints.
    #
    if (filepath):
        with open(filepath) as file:
            condict = {}
            for c in model.getConstrs():
                condict[c.ConstrName] = c

            fcons  = []
            for i, line in enumerate(file):
                tokens = line.split()
                if len(tokens) != 1:
                    print("Warning: bogus line ", i, ", ignoring")
                    continue
                fcons.append(condict[tokens[0]])
            cons  = fcons
    else:
        cons  = mcons

    mcopy     = model.copy()
    constoadd = []
    mconsdict = {}
    for c in mcons:
        mconsdict[c.ConstrName] = c
        
    for c in mcopy.getConstrs():
        if c.ConstrName in mconsdict:
            continue
        constoadd.append((mcopy.getRow(c), c.sense, c.rhs))
        mcopy.remove(c)

    #
    # mcopy now contains the constraints of interest for the first solve.
    #
    mcopy.optimize()
    for c in constoadd:
        mcopy.addConstr(c[0], c[1], c[2])
    mcopy.setParam("Method", 1)      # Make use of dual feasible basis.
    mcopy.optimize()
                     

#
#   Takes a list of decimal values and converts uses the Fraction
#   package to provide a rational approximation.
#
def converttofractions(vals):
    for v in vals:
        if not isinstance(v, float):
            print("Value ", v, " is not a float.  Cannot convert.")
            continue
        frac = Fraction(v).limit_denominator()
        print(f"The approx fraction of {v} is {frac}.")

#
#
#
def fixtostartvalues(model, startfile):
    model.read(startfile)
    vars = m.getVars()
    for v in vars:
        v.LB = v.Start
        v.UB = v.Start
    model.update()

NEGLOWER = 1
POSLOWER = 2
ALLLOWER = 3
#
#   Shifts variables with finite but nonzero lower bounds to 0.
#   Need to adjust the rhs for each bound shift by either adding
#   the positive shift value times the matrix column to the rhs in
#   the case of negative lower bound, and subtracting the positive shift
#   value times the matrix column in the case of a positive lower bound.
#   varstoshift == None (default) means apply shifting to all variables;
#   else input a list of variables to which the shifting is applied.
#
#   Alters the model, so use a copy if you want to preserve the original
#   model.
#
#   TODO: cache the rhs changes yourself so you don't have to call model.update()
#   in each pass through the loop.
#
def shift_lowerbounds(model, bdstoshift, varstoshift=None):
    doneg  = bdstoshift == NEGLOWER or bdstoshift == ALLLOWER
    dopos  = bdstoshift == NEGLOWER or bdstoshift == ALLLOWER
    doboth = bdstoshift == ALLLOWER

    if varstoshift == None:
        vars   = model.getVars()
    else:
        vars = varstoshift
    inf    = math.inf
    objchg = 0.0
    if doboth:
        for v in vars:
            if v.LB > -inf and v.LB < inf and v.LB != 0:
                objchg = shiftonelowerbound(model, v, objchg)
    elif dopos:
        for v in vars:
            if v.LB < inf and v.LB > 0:
                objchg = shiftonelowerbound(model, v, objchg)
    elif doneg:
        for v in vars:
            if v.LB > -inf and v.LB < 0:
                objchg = shiftonelowerbound(model, v, objchg)
    else:
        print("Invalid bound shift ID: ", bdstoshift, "; no bounds shifted")
        return(-1)
    model.addVar(lb=objchg, ub=objchg, obj=1.0, vtype=GRB.CONTINUOUS, \
                 name="ConstantTerm")
    model.update()
                
def shiftonelowerbound(model, v, objchg):
    delta = -v.LB
    col   = model.getCol(v)
    for i in range(col.size()):
        coeff    = col.getCoeff(i)
        con      = col.getConstr(i)
        con.rhs += delta*coeff
    v.LB  = 0.0
    v.UB += delta
    model.update()
    objchg -= delta*v.Obj
    return objchg

#
#   Split free variables into difference of two nonnegative variables.
#   Also works on other variables whose lower bound is the opposite of
#   its upper bound.  varstosplit identifies the candidate list of variables
#   to split.  If not set or None, will treat all variables as candidates.
#   Only those that meet the "mirrored variable" criteria will
#   be split.   Most common use case will be to split free variables.
#
#   Alters the model, so use a copy if you want to preserve the original
#   model.
#
#   Tests:  1) Test on a model with only QCs.
#           2) Test on a model with linear obj terms but no linear constraints.
#
def split_mirroredvars(model, varstosplit=None):
    model.update()
    if varstosplit == None:
        varstosplit = model.getVars()

    newvarlist = []           # This var and dict is used to handle
    newvardict = {}           # quadratic objective and constraints
    for v in varstosplit:
        if v.UB != -v.LB:
            continue
        if v.UB == 0.0 and v.LB == 0.0:       # No need to split vars
            continue                          # fixed at zero
        #
        #  We have a candidate variable to split.
        #
        bnd        = v.UB
        varname    = v.VarName 
        chgvarname = v.VarName  + "_GRBPlus"
        v.VarName  = chgvarname
        v.LB       = 0
        col        = model.getCol(v)
        splitcol   = gp.Column()
        coeflist   = []
        conlist    = []
        for k in range(col.size()):
            coeflist.append(-col.getCoeff(k))
            conlist.append(col.getConstr(k))
        if col.size() > 0:
            splitcol.addTerms(coeflist, conlist)
        newvar = model.addVar(lb=0.0, ub=bnd, obj=-v.obj, vtype=v.Vtype, \
                              name=varname + "_GRBMinus", column=splitcol)
        newvarlist.append(newvar)
        newvardict[chgvarname] = newvar

    #
    # Linear constraints completed; now update any QCs that contain
    # mirrored variables that need to be split.
    #
    for qcon in model.getQConstrs():
        quadexpr = model.getQCRow(qcon)
        split_quadexpr(quadexpr, newvardict)
        #
        # Quadratic expression is updated; need to create a new
        # QC with the update quadexpr and delete the old one.
        #
        model.addQConstr(quadexpr, qcon.QCSense, qcon.QCRHS, qcon.QCName)
        model.remove(qcon)
                    
    #
    # Quadratic constraints completed; now the quadratic objective if it 
    # contains any mirrored variables that need to be split.
    #
    objexpr = model.getObjective()
    if isinstance(objexpr, gp.QuadExpr):
        split_quadexpr(objexpr, newvardict)
        model.setObjective(objexpr)
    model.update()
    return newvardict
#
#   Takes a quadratic expression and replaces mirrored variables with
#   the appropriate difference of two nonnegative variables.  All
#   mirrored variables are provided in newvardict, which maps the
#   original mirrored variable to its added variable.
    
def split_quadexpr(quadexpr, newvardict):
    for k in range(quadexpr.size()):
        xi = quadexpr.getVar1(k)
        xj = quadexpr.getVar2(k)
        qcoef = quadexpr.getCoeff(k)
        if xi.varName in newvardict:      # xi is split into xi(+) - xi(-)
            # xi(-) * xj term
            quadexpr.addTerms(-qcoef, xj, newvardict[xi.VarName])
            if xj.varName in newvardict: # both xi and xj split
                # xi * xj(-) term
                quadexpr.addTerms(-qcoef, xi, newvardict[xj.VarName])
                # xi(-) * xj(-) term
                quadexpr.addTerms(qcoef, newvardict[xi.VarName], \
                                      newvardict[xj.VarName])
        elif xj.varName in newvardict: # xj is split, but xi is not.
            # xi * xj(-) term
            quadexpr.addTerms(-qcoef, xi, newvardict[xj.VarName])


    
#
#   TODO utilities
#   2) single to multi objective. Kostja already has?
#   5) Function moreStats(model) to break down variable and constraint
#      types.   Kostja?


