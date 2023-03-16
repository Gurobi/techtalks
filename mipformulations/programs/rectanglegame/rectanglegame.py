import sys
import copy
import pygame
import math
import gurobipy as gp
from gurobipy import GRB

WHITE    = (255,255,255)
BLACK    = (0, 0, 0)
RED      = (255, 0, 0)
BLUE     = (0, 100, 200)
GREEN    = (0, 255, 0)
GREY     = (128, 128, 128)

WIDTH    = 1200
HEIGHT   = 700
SLACK    = 0
BUTTONS  = 2
RESET    = 0
QUIT     = 1


MINRANDRECTS   = 5
MAXRANDRECTS   = 10
RECTCOLOR      = BLUE
OPTCOLOR       = GREEN
DEVILCOLOR     = GREY
BORDERWIDTH    = 2
BORDERCOLOR    = BLACK
SELECTEDCOLOR  = RED
SELECTEDWIDTH  = 5

CURRENT  = 0
PREVIOUS = 1
GUROBI   = 2
DEVIL    = 3

ERROR    = -1
OK       = 0


#
#   Input rectangles from file, or generate them randomly if
#   no file is specified.   Each line of the file
#   contains a the height and width of a single rectangle
#   These two values are used to initialize the rectangle
#   as an instance of the rectangle class configured above.
#   Also check if we need to shrink the rectangles because
#   they will potentially overflow the screen when concatenated
#   vertically or horizontally.  TODO: add a read_rectangles
#   method to randomly generate a specified number of rectangles.
#
def get_rectangles(rectangles, filename, regionheight, regionwidth, \
                   numrects=None):
    totarea = 0
    if filename == None:   # Generate rectangles randomly
        import random
        if numrects != None:
            numrects = random.randint(MINRANDRECTS, MAXRANDRECTS)
            #
            # Figure out min and maximum rectangle heights and widths based
            # on number of rectangles and game region height and width.
            # Try a simple heuristic approach; sum of rectangle widths and  
            # heights takes up half of the game region width and height
            #
            sumheight = 0
            sumwidth  = 0
            avgheight = int(regionheight/numrects)
            minheight = int(avgheight/8)
            maxheight = 4*minheight
            avgwidth  = int(regionwidth/numrects)
            minwidth  = int(avgwidth/8)
            maxwidth  = 4*minwidth
            filename  = "rectangles_" + str(numrects) + "_" + str(avgheight) + \
                        "_" + str(avgwidth)
            with open(filename + ".txt", "w") as f:
                for j in range(numrects):
                    intheight = random.randint(minheight, maxheight)
                    intwidth  = random.randint(minwidth, maxwidth)
                    sumheight += intheight
                    sumwidth  += intwidth
                    totarea   += float(intwidth*intheight)
                    thisrect  = pygame.Rect(0, 0, intwidth, intheight)
                    rectangles.append(thisrect)
                    f.write(str(intheight) + "  " + str(intwidth) + "\n")
           
    else:    
        with open(filename, "r") as f:
            rects = f.readlines()
            f.close()
        num       = len(rects)
        sumheight = 0
        sumwidth  = 0
        import re
        for j in range(num):
            onerect   = re.split(' |, |; |\n', rects[j])
            onerect   = [x for x in onerect if x]
            if len(onerect) != 2:
                print("Ignoring line ", j, "; unexpected input.")
                print(rects[j])
                continue
            intheight = int(onerect[0])
            intwidth  = int(onerect[1])
            sumheight += intheight
            sumwidth  += intwidth
            totarea   += float(intwidth*intheight)
            thisrect  = pygame.Rect(0, 0, intwidth, intheight)
            rectangles.append(thisrect)
    #
    # If the rectangles won't fit easily on the screen scale them
    # so that they do.   Do likewise if they are too small.  Their
    # heights and widths relative to each other will remain unchanged.
    #
    ratio    = totarea / float(WIDTH*HEIGHT)
    scalefac = 1.0
    if ratio < .05 or ratio > .20:
        #
        # Rescale so that rectangles comprise 10 percent of area of screen.
        #
        scalefac  = math.sqrt(.10/ratio)
        for rect in rectangles:
            rect.width  = round(scalefac*rect.width)
            rect.height = round(scalefac*rect.height)
    print("Sum of rectangle areas = ", round(totarea*scalefac*scalefac))

#
#   Initialize the rectangle positions for the game display.
#   Initial x,y coordinates after reading in all rectangles is just
#   (0,0).  Create a feasibility MILP to scatter than about the
#   game playing region of the screen so that they don't overlap,
#   and at least one rectangle touches each of the 4 boundaries
#   of the region.
#
def init_rectangle_positions(rectangles, gameheight, gamewidth, gamex, gamey):
    nrects = len(rectangles)
    m      = gp.Model("screeninit_" + str(nrects))
    #
    # Variables for model.  (x,y) determine positions of upper left
    # corner of each rectangle.   Binaries Zx and Zy are used to
    # enforce conditions that ensure that at least one rectangle
    # touches each of the 4 boundaries of the game region.
    #
    x   = m.addVars(nrects, ub=gamewidth, vtype=GRB.INTEGER, name="x")
    y   = m.addVars(nrects, ub = gameheight, vtype=GRB.INTEGER, name="y")
    Zx  = m.addVars(nrects, vtype=GRB.BINARY, name="Zx")
    Zy  = m.addVars(nrects, vtype=GRB.BINARY, name="Zy")
    Ux  = m.addVars(nrects, vtype=GRB.BINARY, name="Ux")
    Uy  = m.addVars(nrects, vtype=GRB.BINARY, name="Uy")
    #
    # Constraints that ensure that at least one rectangle touches
    # the boundary of the game region in the initial layout.  This means:
    # 1) At least one rectangle has an x coordinate of 0
    # 2) At least one rectangle has a y coordinate of 0
    # 3) At least one rectangle has x + w = screen width (boundary at right)
    # 4) At least one rectangle has y + h = screen height (boundary at bottom)
    #
    for j in range(nrects):
        width  = rectangles[j].width
        height = rectangles[j].height
        # Constraints 1, 2
        m.addConstr(x[j] - gamewidth * Zx[j] <= 0, name="LeftRect_" + str(j))
        m.addConstr(y[j] - gameheight * Zy[j] <= 0, name="TopRect_" + str(j))
        # Constraint 3
        wbar    = gamewidth - width
        x[j].UB = wbar
        m.addConstr(x[j] - wbar * Ux[j] >= 0, name="RightRect_" + str(j))
        # Constraint 4
        hbar    = gameheight - height
        y[j].UB = hbar
        m.addConstr(y[j] - hbar * Uy[j] >= 0)
    #
    # Constraints involving individual variables for conditions 1-4 above.
    # Now add the cumulative conditions on the binaries Zx, Zy, Ux and Uy
    # that enforce the "at least" conditions for each of the constraints
    #
    m.addConstr(Zx.sum() <= nrects - 1, name="LeftBoundary")
    m.addConstr(Zy.sum() <= nrects - 1, name="TopBoundary")
    m.addConstr(Ux.sum() >= 1, name="RightBoundary")
    m.addConstr(Uy.sum() >= 1, name="BottomBoundary")
    #
    # Constraints to ensure no overlap of any pair of rectangles.  This means:
    # at least one of the following 4 conditions must be true for each pair
    # of rectangles:
    # 5) rectangle i is to the left of rectangle j
    # 6) rectangle i is to the right of rectangle j
    # 7) rectangle i is above rectangle j
    # 8) rectangle i is below rectangle j
    #
    # In addition, to faciliate moving the rectangles around, add some
    # space to the non overlap condition.
    #
    # Additional binaries to enforce no overlap conditions involving
    # width and height, respectively.
    #
    Wxx = m.addVars(nrects, nrects, vtype=GRB.BINARY, name="Wxx")
    Hyy = m.addVars(nrects, nrects, vtype=GRB.BINARY, name="Hyy")
    for i in range(nrects):
        wi = rectangles[i].width
        hi = rectangles[i].height
        for j in range(i+1,nrects):
            wj = rectangles[j].width
            hj = rectangles[j].height
            # Constraints 5, 6, 7 and 8
            m.addConstr(x[j] - x[i] >= wi + SLACK - gamewidth*Wxx[(i,j)], \
                        name="no_overlap_width_" + str(i) + "_" + str(j))
            m.addConstr(x[i] - x[j] >= wj + SLACK - gamewidth*Wxx[(j,i)], \
                        name="no_overlap_width_" + str(j) + "_" + str(i))
            m.addConstr(y[j] - y[i] >= hi + SLACK - gameheight*Hyy[(i,j)], \
                        name="no_overlap_height_" + str(i) + "_" + str(j))
            m.addConstr(y[i] - y[j] >= hj + SLACK - gameheight*Hyy[(j,i)], \
                        name="no_overlap_height_" + str(j) + "_" + str(i))
            #
            # Constraints 5 - 8 are true when associated binary is 0
            # At least one binary must be 0
            #
            m.addConstr(Wxx[i,j] + Wxx[(j,i)] + Hyy[(i,j)] + Hyy[(j,i)] <= 3, \
                        name = "no_overlap_" + str(i) + "_" + str(j))
    m.ModelName = "init_" + str(nrects) + "rects"
    m.write(m.ModelName + ".lp")
    m.optimize()
    #
    # Update the rectangles with their starting positions for the game
    # based on the results of the optimization.  Upon completion, rectangles
    # are ready for display in the game area, which will be done by
    # the calling routine or one of its henchmen.
    #
    if m.status == GRB.OPTIMAL:
        for k in range(nrects):
            rectangles[k].move_ip(gamex + x[k].X, gamey + y[k].X)
        return OK
    else:   # TODO add a dictionary to translate number to message
        print("Unexpected optimization status when initializing: ", m.status)
        return ERROR


#
#   Optimize the rectangle positions for the game display.
#   Given the list of rectangle data, find the positions
#   that minimize the area of the enclosing rectangle without
#   any overlap.
#
def optimize_rectangle_positions(rectangles, gameheight, gamewidth, \
                                 gamex, gamey, scores, relax=False):
    nrects = len(rectangles)
    m      = gp.Model("bestarrangement_" + str(nrects))
    #
    # Variables for model.  (x,y) determine positions of upper left
    # corner of each rectangle.  Wxx and Hyy enforce no overlap conditions
    # involving width and height, respectively.  XStar and YStar measure
    # the width and height of the enclosing rectangle.
    #
    x   = m.addVars(nrects, ub=gamewidth, vtype=GRB.INTEGER, name="x")
    y   = m.addVars(nrects, ub = gameheight, vtype=GRB.INTEGER, name="y")
    if relax:
        Wxx = m.addVars(nrects, nrects, ub=1.0, vtype=GRB.CONTINUOUS, \
                        name="Wxx")
        Hyy = m.addVars(nrects, nrects, ub=1.0, vtype=GRB.CONTINUOUS, \
                        name="Hyy")
    else:    
        Wxx = m.addVars(nrects, nrects, vtype=GRB.BINARY, name="Wxx")
        Hyy = m.addVars(nrects, nrects, vtype=GRB.BINARY, name="Hyy")
    XStar = m.addVar(name="Xstar")
    YStar = m.addVar(name="Ystar")
    #
    # Constraints to ensure no overlap of any pair of rectangles.  This means:
    # at least one of the following 4 conditions must be true for each pair
    # of rectangles:
    # 1) rectangle i is to the left of rectangle j
    # 2) rectangle i is to the right of rectangle j
    # 3) rectangle i is above rectangle j
    # 4) rectangle i is below rectangle j
    #
    for i in range(nrects):
        #
        #  Also specify constraint for the objective variables XStar and
        #  YStar.   Each rectangle defines one constraint for each of those
        #  two variables.
        #
        wi = rectangles[i].width
        hi = rectangles[i].height
        m.addConstr(x[i] + wi <= XStar, name = "DefXStar_" + str(i))
        m.addConstr(y[i] + hi <= YStar, name = "DefYStar_" + str(i))        
        for j in range(i+1,nrects):
            wj = rectangles[j].width
            hj = rectangles[j].height
            # Constraints 1, 2, 3 and 4
            m.addConstr(x[j] - x[i] >= wi - gamewidth*Wxx[(i,j)], \
                        name="no_overlap_width_" + str(i) + "_" + str(j))
            m.addConstr(x[i] - x[j] >= wj - gamewidth*Wxx[(j,i)], \
                        name="no_overlap_width_" + str(j) + "_" + str(i))
            m.addConstr(y[j] - y[i] >= hi - gameheight*Hyy[(i,j)], \
                        name="no_overlap_height_" + str(i) + "_" + str(j))
            m.addConstr(y[i] - y[j] >= hj - gameheight*Hyy[(j,i)], \
                        name="no_overlap_height_" + str(j) + "_" + str(i))
            #
            # Constraints 1 - 4 are true when associated binary is 0
            # At least one binary must be 0
            #
            m.addConstr(Wxx[i,j] + Wxx[(j,i)] + Hyy[(i,j)] + Hyy[(j,i)] <= 3, \
                        name = "no_overlap_" + str(i) + "_" + str(j))
    m.setObjective(XStar*YStar, GRB.MINIMIZE)
    m.ModelName = "bestarrangement_" + str(nrects) + "rects"
    if relax:
        m.ModelName = m.ModelName + "_relaxed"
        
    m.write(m.ModelName + ".lp")
    m.setParam("NonConvex", 2)
    m.optimize()
    #
    # Update the rectangles with their starting positions for the game
    # based on the results of the optimization.  Upon completion, rectangles
    # are ready for display in the game area, which will be done by
    # the calling routine or one of its henchmen.
    #
    if m.SolCount > 0:
        if relax:
            scores[DEVIL] = round(m.ObjVal)
        else:
            scores[GUROBI] = round(m.ObjVal)
        for k in range(nrects):
            rectangles[k].move_ip(gamex + x[k].X, gamey + y[k].X)
        return OK
    else:   # TODO add a dictionary to translate number to message
        print("Unexpected optimization status when initializing: ", m.status)
        return ERROR

    

def draw_rectangles(screen, rectangles, fillcolor, borderwidth=0, \
                    bordercolor=None):
    for thisrect in rectangles:
        if borderwidth > 0:
            draw_one_rectangle(screen, thisrect, bordercolor, borderwidth, \
                               bordercolor)
        pygame.draw.rect(screen, fillcolor, thisrect, width=0)

def draw_one_rectangle(screen, onerect, fillcolor, borderwidth=0, \
                       bordercolor=None):
    if borderwidth > 0:
        pygame.draw.rect(screen, bordercolor, onerect, width=borderwidth)
    elif borderwidth < 0:  # Used to deselect a rectangle
        pygame.draw.rect(screen, WHITE, onerect, width=-borderwidth)
    pygame.draw.rect(screen, fillcolor, onerect, width=0)

#
# Erase one rectangle using the color specified by erasecolor
# Erases border as well if borderwidth < 0; specify negative of
# desired border to specify border thickness.
#
def erase_one_rectangle(screen, onerect, erasecolor, borderwidth=0):
    if borderwidth > 0:
        borderwidth = 0
    draw_one_rectangle(screen, onerect, erasecolor, -borderwidth, erasecolor)


#
#   Draws buttons, which consist of a rectangle with some text inside.
#   Returns the position of the rectangles as a dictionary
#   consisting of the button text and the 4-tuple of the rectangle
#   position, i.e. (x, y, width, height) with x and y the horizontal
#   and vertical coordinates relative to the origin at the upper left
#   corner of the screen.
#    
def draw_buttons(screen, buttoncolors, buttontext, leftmost, rightmost, \
                 bottom):    
    buttondict = {}
    #
    # Find the bottom, left and right most positions of the binary matrices
    #
    vspace     = HEIGHT - (bottom + 1)
    hspace     = rightmost - leftmost
    #
    # Center the buttons within the rectangle with upper left corner at
    # (leftmost+1, bottom+1), width hspace and height vspace
    #
    numbuttons   = len(buttoncolors)
    buttonwidth  = int(hspace/(2*numbuttons))
    buttonheight = int(vspace/3)
    delta        = int((hspace - numbuttons*buttonwidth)/numbuttons)
    x            = leftmost          # first button coordinates
    y            = bottom   + buttonheight
    #
    # Create the position 4-tuples for the rectangles associated with
    # the buttons, then draw them.  + and - start out grey, as they
    # are inactive until a binary matrix is selected, at which point
    # their color will change to blue when they become usable.
    #
    for j in range(numbuttons):
        thisrect  = pygame.Rect(x, y, buttonwidth, buttonheight)
        buttondict[buttontext[j]] = thisrect
        bfontsize = buttonheight - 2
        bfont     = pygame.font.SysFont("menlo", bfontsize)
        pygame.draw.rect(screen, buttoncolors[j], thisrect)
        blabel    = bfont.render(buttontext[j], 1, BLACK)
        screen.blit(blabel, (x,y))
        x += int(3 * buttonwidth / 2) # move right to next button.
        
    return buttondict             # will need this later to change colors


def draw_one_button(screen, buttonrect, buttoncolor, buttontext):
    bfontsize = buttonrect[3] - 2
    bfont     = pygame.font.SysFont("menlo", bfontsize)
    pygame.draw.rect(screen, buttoncolor, buttonrect)
    blabel    = bfont.render(buttontext, 1, BLACK)
    screen.blit(blabel, (buttonrect[0], buttonrect[1]))


#
#   Draws labels and boxes for user score and Gurobi score.
#
def draw_score_labels(screen, xstart, ystart):
    #
    # Calculate spacing.  4 lines of info spaced equally between bottom
    # of game area and bottom of screen.
    #
    y          = ystart + 2
    x          = xstart
    totspace   = HEIGHT - y
    linespace  = int(totspace/6)
    font       = pygame.font.SysFont("menlo", linespace)
    label      = font.render("Scores", 1, BLACK)
    screen.blit(label, (x, y))
    label      = font.render("You :", 1, BLACK)
    y         += int(11 * linespace / 10)
    ylrect     = (x, y, label.get_width(), \
                  label.get_height())    # print scores to the right of this
    screen.blit(label, (x, y))
    surface    = pygame.image.load("Gurobitope.png")
    grbtope    = surface.convert_alpha()
    y         += int(11 * linespace / 10)
    grbrect    = grbtope.get_rect()
    grbrect[0] = x
    grbrect[1] = y
    screen.blit(grbtope, (x, y))
    x         += grbrect[2]
    grblabel   = font.render(":", 1, BLACK)
    w          = grbtope.get_width() + grblabel.get_width()
    h          = grbtope.get_height()
    # GRB score to the right of glrect
    glrect     = (x + grblabel.get_width(), y, w, grblabel.get_height())   
    screen.blit(grblabel, (x,y))
    # Now set up the Devil icon for the configuration from the LP relaxation
    surface    = pygame.image.load("Devil.png")
    devil      = surface.convert_alpha()
    x          = xstart
    y         += 2 * int(11 * linespace / 10)   # a hack.
    dvlrect    = devil.get_rect()
    dvlrect[0] = x
    dvlrect[1] = y
    screen.blit(devil, (x, y))
    x         += dvlrect[2]
    dvllabel   = font.render(" :", 1, BLACK)
    w          = devil.get_width() + dvllabel.get_width()
    h          = devil.get_height()
    # Devil  score to the right of dvlrect
    dvllabelrect = (x + dvllabel.get_width(), y, w, dvllabel.get_height())   
    screen.blit(dvllabel, (x,y))
    #
    # Need to return Gurobitope coordinates since it is clickable.
    # Also return score label rectangles to enable positioning of
    # scoring values.
    #
    return grbrect, dvlrect, ylrect, glrect, dvllabelrect

#
#   Computes the width and height info of the boundary rectangle enclosing
#   a colletion of rectangles.   Specifically, returns the mininum and maximum
#   x and y coordinates on the screen.
#
def boundary_info(rectangles):
    minx     = WIDTH
    miny     = HEIGHT
    maxx     = 0
    maxy     = 0
    
    for rect in rectangles:
        if rect.x < minx:
            minx = rect.x
        if rect.x + rect.width > maxx:
            maxx = rect.x + rect.width
        if rect.y < miny:
            miny = rect.y
        if rect.y + rect.height > maxy:
            maxy = rect.y + rect.height
    return (minx, maxx, miny, maxy)

#
#   Calculate scores and display them in the appropriate location using
#   the rectangles in yourrect and grbrect.  Draw the boundary of the
#   current rectangle placement to illustrate the player's current score
#
def update_scores(screen, scores, yourrect, grbrect, dvlrect, rectangles, \
                  oldlines=None):
    if scores[CURRENT] != None:
        scores[PREVIOUS] = scores[CURRENT]
    newscore = 0
    minx     = WIDTH     # TODO; replace next ~15 lines with boundary_info()
    miny     = HEIGHT
    maxx     = 0
    maxy     = 0
    
    for rect in rectangles:
        if rect.x < minx:
            minx = rect.x
        if rect.x + rect.width > maxx:
            maxx = rect.x + rect.width
        if rect.y < miny:
            miny = rect.y
        if rect.y + rect.height > maxy:
            maxy = rect.y + rect.height
        
    newscore = (maxx - minx)*(maxy - miny)
    #
    # Check if score has degraded compared to previous one.
    #
    scorecolor = None
    scorecolor = RECTCOLOR
    if scores[PREVIOUS] != None:
        if newscore > scores[PREVIOUS]:
            scorecolor = RED
        
    scores[PREVIOUS] = newscore
    #
    # Output to screen.   First clear out old output.
    #
    x         = yourrect[0] + yourrect[2]
    y         = yourrect[1]
    # TODO: fix third argument here, which is a hack.
    pygame.draw.rect(screen, WHITE, [x, y, 2*yourrect[2], yourrect[3]])
    scorefont = pygame.font.SysFont("menlo", yourrect[3])
    output    = scorefont.render(str(newscore), 1, scorecolor)
    screen.blit(output, (x, y))
    # Update Gurobi's score if available. First clear out old output.
    # TODO: clean up use of 2*grbrect[2], which is a hack.
    pygame.draw.rect(screen, WHITE, (grbrect[0], grbrect[1], 2*grbrect[2], \
                                     grbrect[3]))
    if scores[GUROBI] != None:
        scorefont = pygame.font.SysFont("menlo", grbrect[3])
        output    = scorefont.render(str(scores[GUROBI]), 1, OPTCOLOR)
        screen.blit(output, (grbrect[0], grbrect[1]))

    pygame.draw.rect(screen, WHITE, (dvlrect[0], dvlrect[1], 2*dvlrect[2], \
                                     dvlrect[3]))
    # Update Devil's score if available. First clear out old output.
    # TODO: clean up use of 2*grbrect[2], which is a hack.
    if scores[DEVIL] != None:
        scorefont = pygame.font.SysFont("menlo", dvlrect[3])
        output    = scorefont.render(str(scores[DEVIL]), 1, DEVILCOLOR)
        screen.blit(output, (dvlrect[0], dvlrect[1]))

    #
    # Text output updated.  Now illustrate the score by drawing the
    # rectangular boundary around the current configuration of rectangles
    # Erase the previous boundary if it exists.
    #
    if len(oldlines) == 4:
        connect_points(screen, oldlines, WHITE, SELECTEDWIDTH)
        oldlines.clear()
    boundarylines = [((minx, miny), (maxx, miny)), \
                     ((maxx, miny), (maxx, maxy)), \
                     ((maxx, maxy), (minx, maxy)), \
                     ((minx, maxy), (minx, miny))] 
    connect_points(screen, boundarylines, RECTCOLOR, 2)
    oldlines += boundarylines

    
        
        
#
#   Tests if the coordinates in pos are contained by the rectangle,
#   which is represented by the horizonal and vertical positions
#   relative to the origin of (0,0) in the upper left hand corner
#   of the screen, followed by the dimensions of the rectangle.
#
def contained(position, rectangle):
    h      = position[0]
    v      = position[1]
    x      = rectangle.x
    y      = rectangle.y
    width  = rectangle.width
    height = rectangle.height
    if x <= h and h <= x + width and y <= v and v <= y + height:
        result = True
    else:
        result = False
    return result

#
#   Draws the lines in the list, which consists of a sequence of
#   (x,y) coordinates of the points to connect.
#
def connect_points(screen, lines, color, thickness):
    firstpt = lines[0]
    lastpt  = lines[-1]
    for line in lines:
        pygame.draw.line(screen, color, line[0], line[1], thickness)

    

def main():
#    import pdb; pdb.set_trace()   
    filename = None
    numrects = None
    if len(sys.argv) == 2:
        thisarg = sys.argv[1]
        if thisarg.isnumeric():
            numrects = int(thisarg)
        else:
            filename = thisarg
    else:
        filename = "rectangles.txt"
    pygame.init()
    screen  = pygame.display.set_mode((WIDTH, HEIGHT))
    pygame.display.set_caption("GUROBI RECTANGLE PACKING GAME")
    screen.fill(WHITE)
    img       = pygame.image.load('icon.png')
    img       = img.convert_alpha()
    imgrect   = img.get_rect()
    imgheight = img.get_height()
    screen.blit(img, (imgrect[0], imgrect[1]))
    pygame.display.set_icon(img)
    pygame.event.pump()
    pygame.display.flip()
    #
    # Need to figure out location and dimensions of game region before
    # we can draw the rectangles.   Can use full width of screen, but
    # height of region must reduce height of screen by the height of the
    # Gurobi icon at the top, and spaces for 3 lines of scoring info and
    # the game buttons at the bottom.
    #
    surface    = pygame.image.load("Gurobitope.png")
    grbtope    = surface.convert_alpha()
    grbrect    = grbtope.get_rect()
    topeheight = grbrect.height
    gameheight = HEIGHT - (imgheight + 3*topeheight)
    gamewidth  = WIDTH
    gamex      = 0
    gamey      = imgheight
    gamerect   = pygame.Rect(gamex, gamey, gamewidth, gameheight)
    rectangles = []
    data_rectangles = []
    get_rectangles(data_rectangles, filename, gameheight, gamewidth, numrects)
    rectangles = copy.deepcopy(data_rectangles)
    status = init_rectangle_positions(rectangles, gameheight, gamewidth, \
                                      gamex, gamey)
    if status == ERROR:
        font       = pygame.font.SysFont("menlo", int(HEIGHT/20))
        label      = font.render("Error; rectangles don't fit.", 1, BLACK)
        x          = int(WIDTH/4)
        y          = int(HEIGHT/2)
        screen.blit(label, (x, y))
    else:     # Blue rectangle with black border
        draw_rectangles(screen, rectangles, RECTCOLOR, borderwidth=3, \
                        bordercolor=BLACK)
        
    labelx = 5
    labely = gameheight + imgheight + 5
    grbtoperect, devilrect, yourlabelrect, grblabelrect, devillabelrect = \
        draw_score_labels (screen, labelx, labely)

    buttoncolors = [BLUE, BLUE]
    buttontext   = ["reset", "quit "]    
    buttondict   = draw_buttons(screen, buttoncolors, buttontext, \
                                int(WIDTH/2), WIDTH, \
                                HEIGHT - (3*topeheight + 1))
    cutoff       = gameheight*gamewidth
    scores       = [cutoff, cutoff, None, None] # Current, previous, GRB, Devil
    oldlines    = []
    update_scores(screen, scores, yourlabelrect, grblabelrect, devillabelrect,\
                  rectangles, oldlines)
    pygame.display.flip()
#
#   Game board setup completed.  Enter the main game loop.   User
#   can select a rectangle and move it to non overlapping positions
#   with other rectangles
#
#   One additional action consists of clicking on the Gurobi polytope
#   in the scores section.   This will results in the formulation and
#   optimization of the model, with the optimal score printed.   TODO:
#   Enable overlay of optimal rectangle configuration on top of user's
#   configuration.  
#
    origrectangles= rectangles.copy()      # Save in case of reset.
    selected_rect = None
    prevsel_rect  = None
    start_rect    = None
    selected      = False
    dragging      = False
    Quit          = False
    oloopcount    = 0                      # for debugging
    while not Quit:
        event = pygame.event.wait()        # one action at a time processed
        if event.type == pygame.MOUSEBUTTONDOWN:
            #
            # Check if user clicked on something of interest.
            # Start with the 2 action buttons
            #
            pos = pygame.mouse.get_pos()
            button_clicked = False
            for btext in buttontext:
                if contained(pos, buttondict[btext]):
                    if btext == buttontext[RESET]:        # Reset button
                        erase_one_rectangle(screen, gamerect, WHITE, \
                                            -SELECTEDWIDTH)
                        # pygame.draw.rect(screen, WHITE, gamerect, 0)
                        # pygame.display.flip()   # debug only
                        button_clicked = True
                        selected_rect  = None
                        prevsel_rect   = None
                        start_rect     = None
                        selected       = False
                        dragging       = False
                        scores         = [cutoff, cutoff, None, None]  
                        rectangles     = copy.deepcopy(origrectangles)
                        update_scores(screen, scores, yourlabelrect, \
                                      grblabelrect, devillabelrect, rectangles,\
                                      oldlines)
                        draw_rectangles(screen, rectangles, RECTCOLOR, \
                                        borderwidth=3, bordercolor=BLACK)
                        pygame.display.flip()

                    elif btext == buttontext[QUIT]:   # Quit button
                        Quit           = True
                        button_clicked = True
                        break
            if button_clicked:
                continue
            #
            # Check if user clicked on the Gurobitope to request the
            # optimal configuration.
            #
            if contained(pos, grbtoperect):
                optrectangles = copy.deepcopy(data_rectangles)
                status = optimize_rectangle_positions(optrectangles, \
                                                      gameheight, gamewidth, \
                                                      gamex, gamey, scores)
                if status == ERROR:
                    font       = pygame.font.SysFont("menlo", int(HEIGHT/20))
                    errmsg     = "Error; unexpected solve status."
                    label      = font.render(errmsg, 1, BLACK)
                    x          = int(WIDTH/4)
                    y          = int(HEIGHT/2)
                    screen.blit(label, (x, y))
                else:
                    #
                    # Position the optimal configuration as far to
                    # the right of the screen as possible to
                    # keep it separate from the user's configuration
                    #
                    minx, maxx, miny, maxy = boundary_info(optrectangles)
                    width  = maxx - minx
                    shiftx = gamex + gamewidth - maxx
                    for optrect in optrectangles:
                        optrect.x += shiftx 
                    draw_rectangles(screen, optrectangles, OPTCOLOR, \
                                    borderwidth=3, bordercolor=BLACK)
                    #
                    # Draw the boundary rectangle
                    #
                    minx += shiftx
                    maxx += shiftx
                    boundarylines = [((minx, miny), (maxx, miny)), \
                                     ((maxx, miny), (maxx, maxy)), \
                                     ((maxx, maxy), (minx, maxy)), \
                                     ((minx, maxy), (minx, miny))] 
                    connect_points(screen, boundarylines, BLACK, 2)
                    
                    # TODO: clean up use of 2*grbrect[2], which is a hack.
                    pygame.draw.rect(screen, WHITE, (grblabelrect[0], \
                                                     grblabelrect[1], \
                                                     2*grblabelrect[2], \
                                                     grblabelrect[3]))
                    scorefont = pygame.font.SysFont("menlo", grblabelrect[3])
                    output    = scorefont.render(str(scores[GUROBI]), 1, \
                                                 OPTCOLOR)
                    screen.blit(output, (grblabelrect[0], grblabelrect[1]))
                pygame.display.flip()
                
            #
            # Check if user clicked on Dr Evil to request the
            # optimal configuration without the overlap requirements.
            #
            if contained(pos, devilrect):
                devilrectangles = copy.deepcopy(data_rectangles)
                status = optimize_rectangle_positions(devilrectangles, \
                                                      gameheight, gamewidth, \
                                                      gamex, gamey, scores, \
                                                      relax=True)
                if status == ERROR:
                    font       = pygame.font.SysFont("menlo", int(HEIGHT/20))
                    errmsg     = "Error; unexpected solve status."
                    label      = font.render(errmsg, 1, BLACK)
                    x          = int(WIDTH/4)
                    y          = int(HEIGHT/2)
                    screen.blit(label, (x, y))
                else:
                    #
                    # Position the optimal relaxed configuration
                    # in the middle of the screen.  We know that the
                    # Gurobi optimal MIP configuration is on the right, so
                    # just move the player's configuration to the left
                    # of the game screen to make as much space available
                    # as possible.
                    #
                    #
                    # First move player rectangles to the left.
                    #
                    minx, maxx, miny, maxy = boundary_info(rectangles)
                    shift = minx - gamex
                    for rect in rectangles:
                    #    erase_one_rectangle(screen, rect, WHITE, -BORDERWIDTH)
                        rect.x -= shift
                    pygame.draw.rect(screen, WHITE, (minx, miny, maxx - minx, \
                                                     maxy - miny))
                    
                    update_scores(screen, scores, yourlabelrect, grblabelrect, \
                                  devillabelrect, rectangles, oldlines) 
                    draw_rectangles(screen, rectangles, RECTCOLOR, \
                                    borderwidth=3, bordercolor=BLACK)
                    #
                    # Now show the optimal relaxed configuration.
                    #
                    minx, maxx, miny, maxy = boundary_info(devilrectangles)
                    shiftx = (gamewidth - (maxx - minx))/2
                    for dvlrect in devilrectangles:
                        dvlrect.x += shiftx 
                    draw_rectangles(screen, devilrectangles, DEVILCOLOR, \
                                    borderwidth=3, bordercolor=BLACK)
                    minx += shiftx
                    maxx += shiftx
                    #
                    # Create the boundary rectangle for the relaxation
                    # configuration.
                    #
                    boundarylines = [((minx, miny), (maxx, miny)), \
                                     ((maxx, miny), (maxx, maxy)), \
                                     ((maxx, maxy), (minx, maxy)), \
                                     ((minx, maxy), (minx, miny))] 
                    connect_points(screen, boundarylines, BLACK, 2)
                    # TODO: clean up use of 2*grbrect[2], which is a hack.
                    pygame.draw.rect(screen, WHITE, (devillabelrect[0], \
                                                     devillabelrect[1], \
                                                     2*devillabelrect[2], \
                                                     devillabelrect[3]))
                    scorefont = pygame.font.SysFont("menlo", devillabelrect[3])
                    output    = scorefont.render(str(scores[DEVIL]), 1, \
                                                 DEVILCOLOR)
                    screen.blit(output, (devillabelrect[0], devillabelrect[1]))

                    
                pygame.display.flip()
                
            #      
            # No button clicked.  Check if one of the rectangles was selected.
            #
            selected      = False
            selected_rect = None
            for rect in rectangles:
                if contained(pos, rect):
                    selected      = True
                    selected_rect = rect
                    break
                
            if selected:
                if prevsel_rect != None:
                    if selected_rect != prevsel_rect:
                        #
                        # Different rectangle selected; unselect the 
                        # rectangle previously selected (use negative
                        # width argument to deselect).  Newly selected
                        # rectangle will be done subsequently.
                        #
                        erase_one_rectangle(screen, prevsel_rect, \
                                            WHITE, -SELECTEDWIDTH)
                        draw_one_rectangle(screen, prevsel_rect, RECTCOLOR)
                        
                    # Else clicked in already chosen rectangle; no op.
                    # Next operations (button up or dragging started)
                    # handled in subsequent code
                    
                prevsel_rect = selected_rect
                #
                # Highlight the selected rectangle and prepare
                # for moving it.  
                #
                draw_one_rectangle(screen, selected_rect, RECTCOLOR, \
                                   SELECTEDWIDTH, SELECTEDCOLOR)
                pygame.display.flip()
                start_rect = selected_rect
                start_drag = True
                drag_rect  = start_rect
                nextevent  = pygame.event.wait()
                olaprects  = []
                for rect in rectangles:
                    if rect != start_rect:
                        olaprects.append(rect)
                        
                if nextevent.type == pygame.MOUSEBUTTONUP:
                    #
                    # This should only happen in the outer loop when
                    # user selects a rectangle and releases the mouse
                    # botton without dragging, or clicks on an empty area
                    # in the screen.  Go back to main outer loop; wait for
                    # next mouse click to take next action.
                    #
                    continue
                elif nextevent.type == pygame.MOUSEMOTION:
                    dragging = True
                    while True:
                        nextevent = pygame.event.wait()
                        if nextevent.type == pygame.MOUSEMOTION:
                            #
                            # Update dragged rectangle position
                            # First, delete current rectangle
                            # (including border if the initial selected one)
                            #
                            
                            if start_drag == True:
                                w = -SELECTEDWIDTH
                            else:
                                w = 0
                            erase_one_rectangle(screen, drag_rect, WHITE, w)
                            start_drag = False
                            #
                            # Use current mouse position to update rectangle
                            # as it is being dragged.
                            #
                            pos      = pygame.mouse.get_pos()
                            xoff     = pos[0] - drag_rect.centerx
                            yoff     = pos[1] - drag_rect.centery
                            new_rect = drag_rect.move(xoff, yoff)
                            #
                            # Only allow update to display if new rectangle
                            # doesn't collide with any of the other rectangles
                            # and remains within the game playing rectangle
                            #
                            
                            if gamerect.contains(new_rect) and \
                               new_rect.collidelist(olaprects) == -1:
                                # It's all good
                                erase_one_rectangle(screen, drag_rect, WHITE)
                                draw_one_rectangle(screen, new_rect, RECTCOLOR)
                                drag_rect = new_rect
                                pygame.display.flip()
                                prevsel_rect = None
                                continue
                            else:
                                #
                                # Overlap; sound an alarm to disallow and
                                # continue from the last rectangle
                                #
                                print("\a")
                                dragging = False
                                #
                                # Update the list of rectangles to reflect
                                # the new one.
                                #
                                for j in range(len(rectangles)):
                                    if rectangles[j] == start_rect:
                                        rectangles[j] = drag_rect
                                        break
                                update_scores(screen, scores, yourlabelrect, \
                                              grblabelrect, devillabelrect, \
                                              rectangles, oldlines)
                                draw_one_rectangle(screen, drag_rect, RECTCOLOR)
                                pygame.display.flip()
                                break
                                
                        elif nextevent.type == pygame.MOUSEBUTTONUP:
                            #
                            # Use current mouse position to update rectangle
                            # at the moment the mouse button was released.
                            #
                            pos      = pygame.mouse.get_pos()
                            xoff     = pos[0] - drag_rect.centerx
                            yoff     = pos[1] - drag_rect.centery
                            new_rect = drag_rect.move(xoff, yoff)
                            #
                            # Handle the release of the mouse button
                            # after finishing dragging. Make sure final
                            # resting place doesn't overlap any other
                            # rectangle, then place it.   Otherwise,
                            # return the dragged rectangle to its
                            # original position at time of selection.
                            #

                            if new_rect.collidelist(rectangles) == -1:
                                erase_one_rectangle(screen, drag_rect, WHITE)
                                drag_rect = new_rect
                                draw_one_rectangle(screen, new_rect, RECTCOLOR)
                            else:
                                #
                                # Overlap; sound an alarm to disallow and
                                # continue from the last rectangle
                                #
                                print("\a")
                                new_rect = drag_rect
                            #
                            # Update the list of rectangles to reflect
                            # the new one.
                            #
                            for j in range(len(rectangles)):
                                if rectangles[j] == start_rect:
                                    rectangles[j] = new_rect
                                    break
                            update_scores(screen, scores, yourlabelrect, \
                                          grblabelrect, devillabelrect, \
                                          rectangles, oldlines)
                            draw_one_rectangle(screen, new_rect, RECTCOLOR)
                            pygame.display.flip()
                            dragging = False
                            break
                        else:   # shouldn't arrive here; event queue once
                                # MOUSEMOTION is in queue should only stay
                                # at MOUSEMOTION or change to MOUSEBUTTONUP
                            print("Unexpected flow in inner dragging loop.")
                            break
                
#    input("Press any key to continue.")
                
    pygame.quit()




if __name__ == "__main__":
    sys.exit(main())
