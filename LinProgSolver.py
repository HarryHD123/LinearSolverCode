import math
import copy
import time


# This programme is designed to solve linear programing problems using the simplex method

# --- Admin Functions ---

def read_input():
    "Reads Input.txt and formats the problem"

    # Open and read Input
    with open('Input.txt') as Prob:
        Prob = Prob.readlines()

    # Extract Problem Type, Objective Function and Inequalities
    ObjFuncCopy = False
    IneqCopy = False
    EquaCopy = False
    Dual = False

    ObjFunc = []
    Ineq = []
    IneqList = []
    Equa = []
    EquaList = []

    for Line in Prob:

        # Problem Type
        if Line.startswith('    Problem Type'):
            if "1" in Line:
                ProbType = 1
            elif "0" in Line:
                ProbType = -1
            else:
                print("Invalid Problem Type. Please enter '1' for maximisation or '0' for minimisation")
                exit(1)

        # Dual Mode
        if Line.startswith('    Dual'):
            if "1" in Line:
                Dual = True
                ProbType = ProbType *-1
            elif "0" in Line:
                Dual = False
            else:
                print("Invalid Problem Type. Please enter '1' for Dual simplex or '0' for Primal simplex")
                exit(1)

        # Objective Function
        elif "*" in Line:
            ObjFuncCopy = not ObjFuncCopy
        elif ObjFuncCopy:
            Line = Line.strip("\n")
            Line = Line.replace(" ", "")
            ObjFunc = Line

        # Inequalities
        elif "{" in Line:
            IneqCopy = True
        elif "}" in Line:
            IneqCopy = False
        elif IneqCopy:
            Line = Line.strip("\n")
            Line = Line.replace(" ", "")
            IneqList = IneqList + [Line]

        # Equalities
        elif "[" in Line:
            EquaCopy = True
        elif "]" in Line:
            EquaCopy = False
        elif EquaCopy:
            Line = Line.strip("\n")
            Line = Line.replace(" ", "")
            EquaList = EquaList + [Line]

    if ObjFunc == ["\n"]:
        print("Please enter the coefficients of the objective function")
        exit(1)
    elif IneqList == ["\n"] or IneqList == [""] or IneqList == [" "]:
        print("No inequalities detected. \nPlease enter the coefficients of the inequalities.")
        exit(1)

    # Format input from strings to lists and integers
    ObjFunc = ObjFunc.split(",")
    for i in IneqList:
        Ineq.append(i.split(","))
    for i in EquaList:
        Equa.append(i.split(","))

    ObjFunc = [float(i) for i in ObjFunc]

    for i in range(len(Ineq)):
        for j in range(len(Ineq[i])):
            if Ineq[i][j] == '<':
                Ineq[i][j] = '1'
            elif Ineq[i][j] == '>':
                Ineq[i][j] = '-1'

    try:
        Ineq = [[float(j) for j in i] for i in Ineq]
    except ValueError:
        print("No input or a non-numeric input to inequality list detected. "
              "\nPlease input only numbers split with commas except for the sign.")
        exit(1)

    try:
        if not not EquaList[0].strip():
            try:
                Equa = [[int(float(j)) for j in i] for i in Equa]
            except ValueError:
                print("Non-numeric input to equality list detected. \nPlease input only numbers split with commas.")
                exit(1)
    except IndexError:
        pass

    return ProbType, ObjFunc, Ineq, Equa, Dual


def write_file(Filename, Text, Mode):
    "Writes to the a file"

    # This structure allows a string, integer, float, list or list of lists to be written successfully
    with open(Filename, Mode) as f:
        if type(Text) is str or type(Text) is int or type(Text) is float:
            f.write(str(Text) + "\n")
        elif type(Text) is list:
            for i in Text:
                if type(i) is str:
                    f.write(str(i) + "   \t")
                elif type(i) is list:
                    f.write("\n")
                    for j in i:
                        f.write(str("{:.2f}".format(j)) + " \t")  # Tableau always reads as 2.d.p
    f.close()


def create_log(LogDict, Headings, ProbType, Dual):
    "Creates the log.txt and output.txt files"

    Headings = Headings
    HeadingsR = Headings + ["ratio"]

    # Finds the number of pivots logged in LogDict
    NumPivots = 0
    for i in LogDict.keys():
        if "Pivot" in i:
            NumPivots += 1

    # Finds the number of optima logged in LogDict
    NumOptima = 0
    for i in LogDict.keys():
        if "Optimum" in i:
            NumOptima += 1

    # Finds the number of XVariables
    NumXVars = 0
    for i in Headings:
        if "x" in i:
            NumXVars += 1

    # Log.txt initiated
    write_file('log.txt', "This file shows the inner workings of the linear programming solver package.\n"
                          "Each step is recorded and displayed here.\n", 'w')

    # Output.txt initiated
    write_file('output.txt', "This file displays the answers found using linear programming solver package.\n", 'w')

    if ProbType == -1:
        write_file('log.txt', "To solve this minimisation problem the tableau is transposed to find the dual\nproblem"
                              " of the primal problem. The dual problem is then maximised to find the\nsolution to the"
                              " primal problem.\n", 'a')
    elif ProbType == -1 and Dual:
        write_file('log.txt', "To solve this minimisation problem the tableau is transposed and the \ndual simplex"
                              " method is used to solve it.\n", 'a')

    # Logs initial tableau
    write_file('log.txt', "This is the initial Tableau:", 'a')
    write_file('log.txt', Headings, 'a')
    write_file('log.txt', LogDict["Init"], 'a')

    for i in range(1, NumPivots + 1):
        # Logs the pivot
        write_file('log.txt', "\n\nThis is pivot number: {}".format(i), 'a')

        # Logs the pivot column
        if Dual:
            write_file('log.txt',
                       "\nThe pivot column is column: {}".format("y" + str(LogDict["Pivot" + str(i)][2] + 1)),
                       'a')
        else:
            write_file('log.txt',
                       "\nThe pivot column is column: {}".format("x" + str(LogDict["Pivot" + str(i)][2] + 1)),
                       'a')

        # Logs tableau after adding the ratio column
        write_file('log.txt', "\nThis is the tableau with the ratio column added:", 'a')
        write_file('log.txt', HeadingsR, 'a')
        write_file('log.txt', LogDict["Pivot" + str(i)][1], 'a')

        # Logs the pivot row
        write_file('log.txt', "\n\nThe pivot row is row: {}".format(LogDict["Pivot" + str(i)][3] + 1), 'a')

        # Logs tableau after pivot
        write_file('log.txt', "\nThis is the tableau after pivot:", 'a')
        write_file('log.txt', Headings, 'a')
        write_file('log.txt', LogDict["Pivot" + str(i)][0], 'a')

    # Logs answers to both log.txt and output.txt
    write_file('log.txt', "\n\nAn optimum has been found:", 'a')
    write_file('log.txt', "\nThe answers are:", 'a')

    for j in Headings[:-1]:
        try:
            write_file('log.txt', "{} = {:.2f}".format(j, LogDict["Optimum1"][1][j]), 'a')
            write_file('output.txt', "{} = {:.2f}".format(j, LogDict["Optimum1"][1][j]), 'a')
        except KeyError:
            write_file('log.txt', "{} = 0".format(j), 'a')
            write_file('output.txt', "{} = 0".format(j), 'a')

    for i in range(2, NumOptima + 1):  # If there are multiple optima, logs them

        # Multiple optima detection alert
        write_file('log.txt', "\nMultiple optima detected.", 'a')
        write_file('output.txt', "\nMultiple optima detected\n", 'a')

        # Logs the pivot
        write_file('log.txt', "\nThis is multiple optima pivot number: {}".format(i - 1), 'a')

        # Logs the pivot column
        if Dual:
            write_file('log.txt',
                       "\nThe pivot column is column: {}".format("y" + str(LogDict["Pivot" + str(i)][2] + 1)), 'a')
        else:
            write_file('log.txt',
                       "\nThe pivot column is column: {}".format("x" + str(LogDict["Pivot" + str(i)][2] + 1)), 'a')

        # Logs Tableau after adding the ratio column
        write_file('log.txt', "\nThis is the tableau with the ratio column added:", 'a')
        write_file('log.txt', HeadingsR, 'a')
        write_file('log.txt', LogDict["OptPiv" + str(i - 1)][1], 'a')

        # Logs the Pivot Row
        write_file('log.txt', "\n\nThe pivot row is row: {}".format(LogDict["OptPiv" + str(i - 1)][3] + 1), 'a')

        # Logs Tableau after pivot
        write_file('log.txt', "\nThis is the tableau after pivot:", 'a')
        write_file('log.txt', Headings, 'a')
        write_file('log.txt', LogDict["OptPiv" + str(i - 1)][0], 'a')

        # Logs answers to both log.txt and output.txt
        write_file('log.txt', "\n\nAnother optimum has been found:", 'a')
        write_file('log.txt', "\n\nThe answers are:", 'a')
        for j in Headings[:-1]:
            try:
                write_file('log.txt', "{} = {:.2f}".format(j, LogDict["Optimum" + str(i)][1][j]), 'a')
                write_file('output.txt', "{} = {:.2f}".format(j, LogDict["Optimum" + str(i)][1][j]), 'a')
            except KeyError:
                write_file('log.txt', "{} = 0".format(j), 'a')
                write_file('output.txt', "{} = 0".format(j), 'a')

        # Outputs general solution
        l1 = []
        l2 = []
        for j in Headings[:NumXVars]:
            try:
                l1 = l1 + [LogDict["Optimum1"][1][j]]
            except KeyError:
                l1 = l1 + [0]
            try:
                l2 = l2 + [LogDict["Optimum2"][1][j]]
            except KeyError:
                l2 = l2 + [0]
        write_file('output.txt', "\nThe general solution is:\n", 'a')
        write_file('output.txt', ["("] + Headings[:NumXVars] + [")"], 'a')
        write_file('output.txt', "\n\t\t\t\t=", 'a')
        write_file('output.txt', "\t\t\t  Lambda", 'a')
        write_file('output.txt', ["("] + ["\t({})".format(l1)] + ["\t)"], 'a')
        write_file('output.txt', "\n\t\t\t\t+", 'a')
        write_file('output.txt', "\t\t\t (1-Lambda)", 'a')
        write_file('output.txt', ["("] + ["\t({})".format(l2)] + ["\t)"], 'a')
        write_file('output.txt', "\nWhere 0 < Lambda < 1", 'a')

    # Success Message
    if ProbType == 1:
        print("The problem has been solved and maximised successfully\n")
    elif ProbType == -1:
        print("The problem has been solved and minimised successfully\n")


# --- Solver Functions ---

def find_dual(ObjFuncCoeff, IneqCoeff, EquaCoeff):
    "Finds dual of a problem"

    # Dual is found
    TransIneqCoeff = []
    TransObjFuncCoeff = []

    for i in IneqCoeff:
        TransObjFuncCoeff = TransObjFuncCoeff + [i[-2]*-1]
    for j in range(len(ObjFuncCoeff)):
        TempTransIneqCoeff = []
        for i in IneqCoeff:
            TempTransIneqCoeff = TempTransIneqCoeff + [i[j]]
        TempTransIneqCoeff = TempTransIneqCoeff + [ObjFuncCoeff[j]]
        TempTransIneqCoeff = TempTransIneqCoeff + [IneqCoeff[0][-1]]  # Appends the sign
        TransIneqCoeff.append(TempTransIneqCoeff)

    return TransObjFuncCoeff, TransIneqCoeff, EquaCoeff


def simplex_tableau_former(ProbType,ObjFuncCoeff,IneqCoeff,EquaCoeff, Dual):
    # Standard Simplex method tableau

    # Determines whether any equations have been entered
    Equa = True
    for i in EquaCoeff:
        for j in i:
            if not str(j).strip():
                Equa = False

    # Creates a list of variables (x)
    XVarList = []
    for i in range(len(IneqCoeff[0]) - 2):
        if Dual:
            XVarList = XVarList + ["y" + str(i + 1)]
        else:
            XVarList = XVarList + ["x" + str(i + 1)]

    # Creates a list of slack or surplus variables
    SVarList = []
    for i in range(len(IneqCoeff)):
        if Dual:
            SVarList = SVarList + ["x" + str(i + 1)]
        else:
            SVarList = SVarList + ["s" + str(i + 1)]

    # Joins variable lists for use as top line
    AllVar = XVarList + SVarList
    if Dual:
        GenHeadings = ["P", "sol"]
    else:
        GenHeadings = ["z", "sol"]
    Headings = AllVar + GenHeadings

    # Prepares inequality lines
    for i in range(len(IneqCoeff)):
        IneqCoeff[i].pop(len(XVarList)+1)  # removes sign indicator
        IneqCoeff[i] = IneqCoeff[i] + [0] * (len(SVarList + GenHeadings) - 1)  # Adds zeros for s variables
        IneqCoeff[i].append(IneqCoeff[i].pop(len(XVarList)))# Moves the sol column to the end
        IneqCoeff[i][len(XVarList) + i] = 1  # Sets the s variable to 1 for each respective row

    # Prepares equality lines
    for i in range(len(EquaCoeff)):
        EquaCoeff[i] = EquaCoeff[i] + [0] * (len(SVarList + GenHeadings) - 1)  # Adds zeros for s variables
        EquaCoeff[i].append(EquaCoeff[i].pop(len(XVarList)))  # Moves the sol column to the end

    # Prepares objective function line
    if Dual:
        ObjFuncCoeff = [i * ProbType for i in ObjFuncCoeff]  # Sets objective function to zero
    else:
        ObjFuncCoeff = [i * -ProbType for i in ObjFuncCoeff]  # Sets objective function to zero
    ObjFuncLine = ObjFuncCoeff + [0] * (len(SVarList) + len(GenHeadings))  # Adds zeros for s variables
    ObjFuncLine[-2] = 1 * ProbType  # Sets the z column to 1

    # Create the tableau
    Tableau = []
    for Line in IneqCoeff:
        Tableau = Tableau + [Line]
    if Equa:
        for Line in EquaCoeff:
            Tableau = Tableau + [Line]
    Tableau = Tableau + [ObjFuncLine]

    return Tableau, Headings


def formulate_lp_problem(ProbType, ObjFuncCoeff, IneqCoeff, EquaCoeff, Dual):
    "Formulates linear programming problems into tableau format"

    for i in range(len(IneqCoeff)):
        IneqCoeff[i].append(IneqCoeff[i].pop(len(IneqCoeff[0])-2)) # Moves sign to the end for use in program

    if Dual:
        ObjFuncCoeff_Dual, IneqCoeff_Dual, EquaCoeff_Dual = find_dual(ObjFuncCoeff, IneqCoeff, EquaCoeff)
        Tableau, Headings = simplex_tableau_former(ProbType, ObjFuncCoeff_Dual, IneqCoeff_Dual,
                                                             EquaCoeff, Dual)
    else:
        Tableau, Headings = simplex_tableau_former(ProbType,ObjFuncCoeff,IneqCoeff,EquaCoeff, Dual)

    return Tableau, Headings


def pivot(Tableau, PivColumn, PivRow):
    "Performs a pivot"

    for i in range(len(Tableau)):
        if Tableau[i] != Tableau[PivRow] and Tableau[i][PivColumn] != 0:  # Except for pivot row

            try:
                Multiple = Tableau[i][PivColumn] / Tableau[PivRow][PivColumn]  # Multiple to set value to zero found
            except ZeroDivisionError:
                print(
                    "Zero division error for pivot value detected.\n"
                    "This problem cannot be solved, please check the inputs are correct.\n"
                    "(This error can be often be fixed by making sure the problem type and dual are opposite)).")
                exit(2)
            ManipRow = []  # Manipulated row reset
            for j in range(len(Tableau[i])):
                ManipRow = ManipRow + [(Tableau[i][j] - (Multiple * Tableau[PivRow][j]))]  # Pivot performed
            Tableau[i] = ManipRow  # Tableau updated

    return Tableau


def read_answer(Tableau, Headings, Dual):
    "Reads off answers from a complete tableau"

    AnswerDict = {}

    if Dual:
        XVarCol = [i for i, j in enumerate(Headings) if 'x' in j]
        for i in XVarCol:
            AnswerDict["{}".format(Headings[i])] = Tableau[-1][i]
        AnswerDict["P"] = Tableau[-1][-1]

    else:
        for i in range(len(Tableau[0]) - 1):  # For each column
            NonZeroCount = 0
            TempOneRow = 0

            for j in range(len(Tableau)):  # Checks number of non-zeros in column
                if float(Tableau[j][i]) != 0.0:
                    NonZeroCount += 1
                    TempOneRow = j

            if NonZeroCount == 1:  # If there is only one non-zero, reads off answer
                TempOnePos = [i, TempOneRow]
                TempUVectMulti = Tableau[TempOnePos[1]][TempOnePos[0]] / 1
                for k in range(len(Tableau[TempOnePos[1]])):
                    Tableau[TempOnePos[1]][k] /= TempUVectMulti
                AnswerDict["{}".format(Headings[TempOnePos[0]])] = Tableau[TempOnePos[1]][-1] / Tableau[TempOnePos[1]][
                    TempOnePos[0]]

    return AnswerDict, Tableau


def scale_unit_vector(Tableau):
    "Scales unit vectors"
    for i in range(len(Tableau[0]) - 2):  # For each column except z/P
        NonZeroCount = 0
        TempOneRow = 0

        for j in range(len(Tableau)):  # Checks number of non-zeros in column
            if float(Tableau[j][i]) != 0.0:
                NonZeroCount += 1
                TempOneRow = j

        if NonZeroCount == 1:  # If there is only one non-zero, reads off answer
            TempOnePos = [i, TempOneRow]
            TempUVectMulti = Tableau[TempOnePos[1]][TempOnePos[0]] / 1
            for k in range(len(Tableau[TempOnePos[1]])):
                Tableau[TempOnePos[1]][k] /= TempUVectMulti

    return Tableau


def multiple_optima_finder(Tableau, Headings):
    "Checks for multiple optima"

    MultipleOptimaList = []

    XVarCol = [i for i, j in enumerate(Headings) if 'x' in j]

    # for i in range(XCount):  # For each X column
    for i in XVarCol:  # For each X column
        Column = []
        NonZeroCount = 0
        for j in range(len(Tableau)):
            Column.append(Tableau[j][i])
            if float(Tableau[j][i]) != 0.0:
                NonZeroCount += 1

        if Tableau[-1][i] == 0 and NonZeroCount > 1:
            MultipleOptimaList.append(i)

    return MultipleOptimaList


def lin_prog_solver(Tableau, Headings, ProbType, Dual):
    "Solves linear programming problems"

    # Sets variables to be used in the code
    TableauInit = copy.deepcopy(Tableau)
    LogDict = {}
    LogDict["Init"] = TableauInit

    if Dual:
        ProbType = ProbType * -1

    # Loop variables set
    Solved = False
    PivNum = 0
    Timeout = time.time() + 10

    while not Solved:  # Loops code until a solution is found
        # Timeout added to prevent waiting if there are any errors
        if time.time() > Timeout:
            print("Timeout reached. No optimal solution could not be found.")
            exit(3)

        # Resets loop variables
        PivNum += 1

        # Checks if problem is solved
        if ((ProbType == 1 and not Dual) or (ProbType == -1 and Dual)) and all(i >= 0 for i in Tableau[-1]):
            Solved = True
            break
        elif ((ProbType == -1 and not Dual) or (ProbType == 1 and Dual)) and all(i <= 0 for i in Tableau[-1]):
            Solved = True
            break

        # Find pivot column
        if ProbType == 1 or (ProbType == -1 and Dual):
            PivCol = Tableau[-1].index(min(Tableau[-1]))
        elif (ProbType == -1 and not Dual):
            PivCol = Tableau[-1].index(max(Tableau[-1][:-2]))

        # Add ratios
        Ratios = []
        TableauRatio = copy.deepcopy(Tableau)
        for Row in TableauRatio[:-1]:
            try:
                Ratio = Row[-1] / Row[PivCol]
            except ZeroDivisionError:
                Ratio = math.inf
            except TypeError:
                Ratio = float(Row[-1]) / float(Row[PivCol])
            if Row[PivCol] >= 0:
                Ratios.append(Ratio) # Append used instead of extend to counteract type errors
            else:
                Ratios.append(math.inf) # This fixes the issue of pivoting at -0
            Row.append(Ratio)

        # Find pivot row
        PivRow = Ratios.index(min(i for i in Ratios if i >= 0))

        # Perform pivot
        Tableau = pivot(Tableau, PivCol, PivRow)

        # Scale any unit vectors
        Tableau = scale_unit_vector(Tableau)

        # Log pivot
        LogDict["Pivot" + str(PivNum)] = [copy.deepcopy(Tableau), copy.deepcopy(TableauRatio), PivCol, PivRow]

    # Read off answers
    AnswerDict, Tableau = read_answer(Tableau, Headings, Dual)

    # Log final tableau and answers
    LogDict["Optimum1"] = [copy.deepcopy(Tableau), AnswerDict]

    # Checks for multiple optima
    MultipleOptimaIndexList = multiple_optima_finder(Tableau, Headings)
    OptimaNum = 1

    # For each additional optima found a pivot is needed
    # Log other optima, tableaus and answers
    for i in MultipleOptimaIndexList:
        OptimaNum += 1
        # Add ratios
        Ratios = []
        TableauRatio = copy.deepcopy(Tableau)
        for Row in TableauRatio[:-1]:
            try:
                Ratio = Row[-1] / Row[i]
            except ZeroDivisionError:
                Ratio = math.inf
            except TypeError:
                Ratio = float(Row[-1]) / float(Row[i])
            Ratios.append(Ratio)  # Append used instead of extend to counteract type errors
            Row.append(Ratio)

        # Find pivot row
        PivRow = Ratios.index(min(i for i in Ratios if i >= 0))

        # Perform pivot
        Tableau = pivot(Tableau, i, PivRow)

        # Log pivot
        LogDict["OptPiv" + str(OptimaNum - 1)] = [copy.deepcopy(Tableau), copy.deepcopy(TableauRatio), i, PivRow]

        # Read off answers
        AnswerDict, Tableau = read_answer(Tableau, Headings)

        # Log tableau and answers
        LogDict["Optimum" + str(OptimaNum)] = [copy.deepcopy(Tableau), AnswerDict]

    return LogDict


# --- MAIN ---

def main():
    "The Main function that runs the linear programming solver"
    ProbType, ObjFuncCoeff, IneqCoeff, EquaCoeff, Dual = read_input()  # Input Read
    Tableau, Headings = formulate_lp_problem(ProbType, ObjFuncCoeff, IneqCoeff, EquaCoeff, Dual)  # Tableau Formulated
    LogDict = lin_prog_solver(Tableau, Headings, ProbType, Dual)  # Problem Solved
    create_log(LogDict, Headings, ProbType, Dual)  # Problem Logged


if __name__ == '__main__':
    main()









    '''
LEGACY (PLEASE IGNORE)
    
        Swap sign
        if IneqCoeff[i][-1] == 1 and ProbType == -1:  # If minimisation
            IneqCoeff[i] = [j * -1 for j in IneqCoeff[i]]  # Turns <= into >= so
        elif IneqCoeff[i][-1] == -1 and ProbType == 1:  # If maximisation
            IneqCoeff[i] = [j * -1 for j in IneqCoeff[i]]  # Turns >= into <=
        
    
Read Answer
        if Dual:
            XVarCol = [i for i, j in enumerate(Headings) if 'x' in j]
            for i in XVarCol:
                AnswerDict["{}".format(Headings[i])] = Tableau[-1][i]
            AnswerDict["P"] = Tableau[-1][-1]
    '''