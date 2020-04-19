from utility import *
import txtparser
import model
import pmodel
import re
import copy as cp
from timeit import default_timer as timer


#rename all TF to TFe ID in transcription if TF has effector
# def setTFeffectorID(mod):
#     tfes = {}
#     for tfc in mod.tfconformations.keys():
#         tfe = mod.tfconformations[tfc].tfe
#         tf = mod.tfconformations[tfc].tf
#         tfes[tf] = tfe
#     for tr in mod.transcriptions.keys():
#         for tf in mod.transcriptions[tr].tfs:
#             if tf in tfes:
#                 mod.transcriptions[tr].tfs.remove(tf)
#                 mod.transcriptions[tr].tfs.append(tfes[tf])


#check if gene is linked (direct/indirect) to metnet rxns
# def checkGeneMnet(mod, chked, tgs):
#     for tg in tgs:
#         if tg not in chked:
#             chked[tg] = 1
#             if mod.genes[tg].mnet==True:
#                 return True
#                 break
#             elif mod.genes[tg].tgs:
#                 checkGeneMnet(mod, chked, mod.genes[tg].tgs)
#             else: return False


#check if gene is linked (direct/indirect) to metnet rxns
def checkGeneMnet(mod):
    chkin = mod.genes.keys()
    while chkin:
        rms = []
        for id in chkin:
            if mod.genes[id].mnet:
                for tf in mod.genes[id].tfs:
                    mod.genes[tf].mnet = True
                rms += [id]
        if not rms: break
        else:
            for x in rms: chkin.remove(x)


def readModel(modelfile, id='', name='', ext='_b'):
    f = open(modelfile, 'r')
    mod = model.Model(id=id, name=name)
    try:
        txtparser.parserModelFile(f, mod, ext)
        checkGeneMnet(mod)
        #if mod.tfconformations: setTFeffectorID(mod)
        return mod
    except ParserError, e:
        print 'model syntax error at line ' + e.nline + ': ' + e.msg
    finally:
        f.close()



def readComponents(compfile, mod):
    f = open(compfile, 'r')
    try:
        txtparser.parserCompFile(f, mod)
        #return mod
    except ParserError, e:
        print 'model syntax error at line ' + e.nline + ': ' + e.msg
    finally:
        f.close()


def readExprData(exprfile):
    f = open(exprfile, 'r')
    expr = model.ExpressionData()
    try:
        txtparser.parserExprFile(f, expr)
        return expr
    except ParserError, e:
        print ' problem syntax error at line ' + e.nline + ': ' + e.msg
    finally:
        f.close()


def readProblem(probfile):
    f = open(probfile, 'r')
    prob = model.Problem()
    try:
        txtparser.parserProbFile(f, prob)
        return prob
    except ParserError, e:
        print ' problem syntax error at line ' + e.nline + ': ' + e.msg
    finally:
        f.close()


def writeProblem(self, prob, problemfile):
    pass


class Analysis:
    #minimize' (1) or 'maximize' (-1)
    def __init__(self, mod, obj='', sense=-1, eps=1e-6, options={}, prob='', cmt=''):
        self.mod = mod
        self.analysis = 'fba'
        self.sense = sense
        self.setpmodvar = {} # varid: {lb/ub/cmt: value}, can only set var bounds
        self.setpmodcst = {} # cstid: {expr/lb/ub/cmt: value}, can only add new cst
        self.bookvar = {} #used to book var bounds for renewing
        self.options = options #solver: {option id: value}
        self.cmt = cmt
        self.eps = eps
        self.createPmodel(self.mod)
        if obj: self.setObjective(obj)
        self.prob = prob


    def createPmodel(self, mod):
        self.pmod = pmodel.Pmodel(mod)


    #obj: the string or dict var->coef
    def setObjective(self, obj, sense=-1):
        if obj.__class__ == str: #string
            self.obj = obj
            obj = txtparser.parserLinearExpr(self.obj)
        else: #dict var->coef
            tt = []
            for k, v in obj.items():
                if v == 1:
                    tt += [k]
                else:
                    tt += [str(v) + ' ' + k]
            self.obj = ' + '.join(tt)
        self.sense = sense
        self.pmod.addObjective(obj, self.sense, self.mod)


    #set a variable bounds and book current settings
    def setVarBounds(self, id, lb=float('-Inf'), ub=float('Inf'), cmt=''):
        self.pmod.setVarBounds(id, lb, ub)
        if id in self.setpmodvar:
            self.setpmodvar[id]['lb'] = lb
            self.setpmodvar[id]['ub'] = ub
            if cmt: self.setpmodvar[id]['cmt'] = cmt
        else:
            self.setpmodvar[id] = {'lb':lb, 'ub':ub, 'cmt': cmt}


    #set a expr constraint and book current settings
    def addConstraint(self, id, expr, lb=float('-Inf'), ub=float('Inf'), cmt=''):
        expr2 = txtparser.parserLinearExpr(expr)
        self.pmod.addConstraint(id, expr2, lb, ub)
        if id in self.setpmodcst:
            self.setpmodcst[id]['expr'] = expr
            self.setpmodcst[id]['lb'] = lb
            self.setpmodcst[id]['ub'] = ub
            if cmt: self.setpmodcst[id]['cmt'] = cmt
        else:
            self.setpmodcst[id] = {'expr':expr, 'lb':lb, 'ub':ub, 'cmt': cmt}


    #set a set of variable bounds
    def setVarsBounds(self, csts):
        for id, bds in csts.items():
            self.pmod.setVarBounds(id, bds[0], bds[1])


    #setback to original var bounds and set a new a set of var bounds
    def renewVarsBounds(self, csts):
        if self.bookvar: self.setVarsBounds(self.bookvar)  # set back to initial bounds
        self.bookvar = {} #clear current book
        for id, bds in csts.items():  # set new bounds
            self.bookvar[id] = self.pmod.getVarBounds(id)
            self.pmod.setVarBounds(id, bds[0], bds[1])


    def clearSetpmodvar(self):
        for id in self.setpmodvar:
            if id in self.pmod.vars:
                lb, ub = self.mod.reactions[id].getBounds()
                self.pmod.setVarBounds(id, lb, ub)
            elif id in self.pmod.trvs:
                self.pmod.setVarBounds(id, 0, 1)
            self.setpmodvar.pop(id)


    def clearSetpmodcst(self):
        for id in self.setpmodcst:
            self.pmod.del_component(self.pmod.csts[id])
            self.setpmodcst.pop(id)


    def setSolverOptions(self, solver, option, value=''): #None for default value
        if solver in self.options:
            self.options[solver][option] = value
        else:
            self.options[solver] = {}
            self.options[solver].update({option:value})


    #options="option1=value option2=value" will extend and overwrite the same options in self.options
    # def optsolve(self, solver='glpk', options=None, *args, **kwargs):
    #     self.solver = solver
    #     if options:
    #         for option in options.split():
    #             tt = option.split('=')
    #             self.options[tt[0]] = tt[1]
    #     self.pout = self.pmod.solve(solver, self.options, *args, **kwargs)

    #only use self.options
    def optsolve(self, solver='glpk', *args, **kwargs):
        self.solver = solver
        self.pout = self.pmod.solve(solver, options=self.options.get(solver, None), *args, **kwargs)


    #solver with combination of options
    def optsolveComb(self, solver='glpk2', *args, **kwargs):
        solver0 = re.sub('\d$', '', solver)
        if solver=='glpk2':
            #self.setSolverOptions('glpk', 'tmlim', value=tmlim)
            self.setSolverOptions('glpk', 'dual', value='')
            self.setSolverOptions('glpk', 'nopresol', value='')
            self.optsolve(solver=solver0, *args, **kwargs)
            st = str(self.pout.solver.termination_condition)
            if st != 'optimal':
                self.options['glpk'].pop('nopresol')
                self.optsolve(solver=solver0, *args, **kwargs)
                st = str(self.pout.solver.termination_condition)
                if st != 'optimal':
                    self.options['glpk'].pop('dual')
                    self.optsolve(solver=solver0, *args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    if st != 'optimal':
                        self.setSolverOptions('glpk', 'nopresol', value='')
                        self.optsolve(solver=solver0, *args, **kwargs)
        elif solver=='gurobi2':
            self.setSolverOptions('gurobi', 'FeasibilityTol', value='1e-6')
            self.setSolverOptions('gurobi', 'ScaleFlag', value='-1')
            self.optsolve(solver=solver0, *args, **kwargs)
            st = str(self.pout.solver.termination_condition)
            if st != 'optimal':
                self.setSolverOptions('gurobi', 'ScaleFlag', value='1')
                self.optsolve(solver=solver0, *args, **kwargs)
                st = str(self.pout.solver.termination_condition)
                if st != 'optimal':
                    self.setSolverOptions('gurobi', 'ScaleFlag', value='0')
                    self.optsolve(solver=solver0, *args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    if st != 'optimal':
                        self.setSolverOptions('gurobi', 'FeasibilityTol', value='1e-4')
                        self.setSolverOptions('gurobi', 'ScaleFlag', value='0')
                        self.optsolve(solver=solver0, *args, **kwargs)
        elif solver=='gurobi3':
            self.setSolverOptions('gurobi', 'FeasibilityTol', value='1e-6')
            self.setSolverOptions('gurobi', 'ScaleFlag', value='0')
            self.optsolve(solver=solver0, *args, **kwargs)
            st = str(self.pout.solver.termination_condition)
            if st != 'optimal':
                self.setSolverOptions('gurobi', 'FeasibilityTol', value='1e-4')
                self.setSolverOptions('gurobi', 'ScaleFlag', value='0')
                self.optsolve(solver=solver0, *args, **kwargs)
        else:
            self.optsolve(solver=solver, *args, **kwargs)


    def deactiveConstraints(self):
        pass


    def writePmodel(self, outfile):
        with open(outfile, 'w') as f:
            self.pmod.pprint(ostream=f)


    def writePmodel2(self, outfile):
        self.pmod.writeModel(outfile, self.mod)


    def mprint(self, id):
        self.pmod.myprint(id)

    #show component info.
    def show(self, id):
        pass



class Fba(Analysis):
    def __init__(self, mod, lst='', *args, **kwargs):
        Analysis.__init__(self, mod, *args, **kwargs)
        self.lst = lst
        if lst.__class__ == str: self.lst = lst.split()
        self.results = {} #rxn: flux value
        self.resultsMedia = {}  # results for each medium
        if self.prob: self.setProblem()
        if not self.lst:
            self.lst = self.pmod.rids + self.pmod.gids


    def setProblem(self):
        if self.prob.obj: self.setObjective(self.prob.obj)
        if self.prob.constraints: self.setVarsBounds(self.prob.constraints)
        if not self.lst: self.lst = self.prob.rxns + self.prob.genes


    def solve(self, solver='glpk2', obj='', med='', *args, **kwargs):
        time_start = timer()
        if obj: self.setObjective(obj)
        self.optsolveComb(solver=solver, *args, **kwargs)
        self.results['solver_status'] = str(self.pout.solver.status)
        self.results['solver_termination'] = str(self.pout.solver.termination_condition)
        self.results['obj'] = self.obj
        if self.results['solver_termination'] == 'optimal':
            self.results['objvalue'] = self.pmod.getObjValue()
            self.results['table'] = self.pmod.getVarValues()
            if med:
                print med + '  ' + str(self.results['objvalue'])+'  '+self.results['solver_termination']
            else: print str(self.results['objvalue'])+'  '+self.results['solver_termination']
        else:
            if med:
                print med + '  ' + 'Optimal_Not_Found' + '  ' + self.results['solver_termination']
            else: print 'Optimal_Not_Found' + '  ' + self.results['solver_termination']
        time_end = timer()
        self.results['runtime'] = time_end - time_start


    def solveMedia(self, solver='glpk2', obj='', *args, **kwargs):
        time_start = timer()
        if self.prob and self.prob.media:
            outmd = {}  # multi-medium results
            for md in self.prob.media:
                self.renewVarsBounds(self.prob.media[md])
                self.solve(solver=solver, obj=obj, med=md, *args, **kwargs)
                outmd[md] = cp.deepcopy(self.results)
            self.resultsMedia['outmd'] = outmd
            time_end = timer()
            self.resultsMedia['runtime'] = time_end - time_start
        else: self.solve(solver=solver, obj=obj, *args, **kwargs)


    def write_fba(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: FBA,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.obj])+'\n')
        f.write('#TimeSpent: ' + str(self.results['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o+'='+str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'Value', 'Comment']) + '\n')
        table = self.results['table']
        if self.lst:
            for id in self.lst:
                f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(table[id]), self.mod.getRxnComment(id)])+'\n')
        else:
            for id, v in table.iteritems():
                f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(v), self.mod.getRxnComment(id)])+'\n')
        f.close()


    def writeMedia_fba(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: FBA,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            f.write('-' * 100 + '\n')
            f.write('#Medium: ' + md + '\n')
            f.write('#Objective: ' + ', '.join([str(results['objvalue']), results['solver_termination'], results['obj']]) + '\n')
            f.write('#TimeSpent: ' + str(results['runtime']) + ' (s)' + '\n')
            f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'Value', 'Comment']) + '\n')
            if results['solver_termination'] == 'optimal':
                table = results['table']
                if self.lst:
                    for id in self.lst:
                        f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(table[id]), self.mod.getRxnComment(id)]) + '\n')
                else:
                    for id, v in table.iteritems():
                        f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(v), self.mod.getRxnComment(id)]) + '\n')
        f.close()


    def write_obj(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: OBJ,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.obj]) + '\n')
        f.write('#TimeSpent: ' + str(self.pout.solver.time) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        f.close()


    def writeMedia_obj(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: OBJ,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        f.write('#Columns: ' + ', '.join(['medium', 'objective', 'objvalue', 'solver_status']) + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            if results['solver_termination'] == 'optimal':
                f.write('\t'.join([md, results['obj'], str(results['objvalue']), results['solver_termination']]) + '\n')
            else: f.write('\t'.join([md, results['obj'], 'NA', results['solver_termination']]) + '\n')
        f.close()



class Fva(Analysis):
    def __init__(self, mod, grp='gene', lst='', *args, **kwargs):
        Analysis.__init__(self, mod, *args, **kwargs)
        self.grp = grp
        self.lst = lst
        if lst.__class__ == str: self.lst = lst.split()
        self.results = {}  # single results
        self.resultsMedia = {}  # results for each medium
        if self.prob: self.setProblem()
        if not self.lst:
            if self.grp == 'reaction': self.lst = self.pmod.rids
            elif self.grp == 'gene': self.lst = self.pmod.gids



    def setProblem(self):
        if self.prob.obj: self.setObjective(self.prob.obj)
        if self.prob.constraints: self.setVarsBounds(self.prob.constraints)
        if not self.lst:
            if self.grp == 'reaction':
                if self.prob.rxns:
                    self.lst = self.prob.rxns
                else: self.lst = self.pmod.rids
            elif self.grp == 'gene':
                if self.prob.genes:
                    self.lst = self.prob.genes
                else: self.lst = self.pmod.gids


    # fix obj_var with max-eps, max+eps
    def solve(self, solver='glpk2', obj='', *args, **kwargs):
        time_start = timer()
        if obj: self.setObjective(obj)
        self.optsolveComb(solver=solver, *args, **kwargs)
        self.results['solver_status'] = str(self.pout.solver.status)
        self.results['solver_termination'] = str(self.pout.solver.termination_condition)
        self.results['obj'] = self.obj
        table = {}
        if self.results['solver_termination'] == 'optimal':
            self.results['objvalue'] = self.pmod.getObjValue()
            self.pmod.addConstraintObj(self.eps)
            print str(self.results['objvalue'])+'  '+self.results['solver_termination']
            n = 0
            for id in self.lst:
                if self.pmod.getVarActive(id) == True:
                    self.setObjective(id, sense=1)
                    self.optsolveComb(solver=solver, *args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    table[id] = {'stmin': st}
                    table[id].update({'min': self.pmod.getVarValue(id)})
                    self.setObjective(id, sense=-1)
                    self.optsolveComb(solver=solver, *args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    table[id].update({'stmax': st})
                    table[id].update({'max': self.pmod.getVarValue(id)})
                    n += 1
                    print str(n), id, 'min='+str(table[id]['min'])+':'+table[id]['stmin'], 'max='+str(table[id]['max'])+':'+table[id]['stmax']  # tst
            time_end = timer()
            self.results['table'] = table
            self.results['runtime'] = time_end - time_start
        else:
            print 'Optimal_Not_Found' + '  ' + self.results['solver_termination']


    def solveMedia(self, solver='glpk2', obj='', *args, **kwargs):
        time_start = timer()
        if self.prob and self.prob.media:
            outmd = {}  # multi-medium results
            for md in self.prob.media:
                print '>Medium: '+md
                self.renewVarsBounds(self.prob.media[md])
                self.solve(solver=solver, obj=obj, *args, **kwargs)
                outmd[md] = cp.deepcopy(self.results)
            self.resultsMedia['outmd'] = outmd
            time_end = timer()
            self.resultsMedia['runtime'] = time_end - time_start
        else: self.solve(solver=solver, obj=obj, *args, **kwargs)


    def write(self, resultsfile):
        zcut = 1e-8  # zero cutoff for gene flux
        f = open(resultsfile, 'w')
        f.write('#Analysis: FVA,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.results['obj']]) + '\n')
        f.write('#TimeSpent: ' + str(self.results['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'Comment']) + '\n')
        table = self.results['table']
        for id in table.keys():
            if self.pmod.getVarActive(id) == True:
                activity = '0'
                stmin = table[id]['stmin']
                stmax = table[id]['stmax']
                min = table[id]['min']
                max = table[id]['max']
                if stmin == 'optimal' and stmax == 'optimal':
                    if (min > zcut and max > zcut) or (min < -zcut and max < -zcut):
                        activity = '1'
                    elif abs(min) < zcut and abs(max) <= zcut:
                        activity = '-1'
                f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(min), stmin, str(max), stmax, activity, self.mod.getRxnComment(id)]) + '\n')
        f.close()


    def writeMedia(self, resultsfile):
        zcut = 1e-8  # zero cutoff for gene flux
        f = open(resultsfile, 'w')
        f.write('#Analysis: FVA,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            f.write('-' * 100 + '\n')
            f.write('#Medium: ' + md + '\n')
            f.write('#Objective: ' + ', '.join([str(results['objvalue']), results['solver_termination'], results['obj']]) + '\n')
            f.write('#TimeSpent: ' + str(results['runtime']) + ' (s)' + '\n')
            f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'Comment']) + '\n')
            if results['solver_termination'] == 'optimal':
                table = results['table']
                for id in table.keys():
                    if self.pmod.getVarActive(id) == True:
                        activity = '0'
                        stmin = table[id]['stmin']
                        stmax = table[id]['stmax']
                        min = table[id]['min']
                        max = table[id]['max']
                        if stmin == 'optimal' and stmax == 'optimal':
                            if (min > zcut and max > zcut) or (min < -zcut and max < -zcut):
                                activity = '1'
                            elif abs(min) < zcut and abs(max) <= zcut:
                                activity = '-1'
                        f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(min), stmin, str(max), stmax, activity, self.mod.getRxnComment(id)]) + '\n')
        f.close()



class Ko(Analysis):
    #grp=reaction/gene
    def __init__(self, mod, grp='gene', lst='', *args, **kwargs):
        Analysis.__init__(self, mod, *args, **kwargs)
        self.grp = grp
        self.lst = lst
        if lst.__class__ == str: self.lst = lst.split()
        self.results = {}
        self.resultsMedia = {} #results for each medium
        if self.prob: self.setProblem()
        if not self.lst:
            if self.grp == 'reaction': self.lst = self.pmod.rids
            elif self.grp == 'gene': self.lst = self.pmod.gids


    def setProblem(self):
        if self.prob.obj: self.setObjective(self.prob.obj)
        if self.prob.constraints: self.setVarsBounds(self.prob.constraints)
        if not self.lst:
            if self.grp == 'reaction':
                if self.prob.rxns:
                    self.lst = self.prob.rxns
                else: self.lst = self.pmod.rids
            elif self.grp == 'gene':
                if self.prob.genes:
                    self.lst = self.prob.genes
                else: self.lst = self.pmod.gids


    def solve(self, solver='glpk2', obj='', *args, **kwargs):
        time_start = timer()
        if obj: self.setObjective(obj)
        self.optsolveComb(solver=solver, *args, **kwargs)
        stbm = str(self.pout.solver.termination_condition)
        if stbm == 'optimal':
            bmwt = self.pmod.getObjValue()
            self.results['solver_termination'] = stbm
            self.results['obj'] = self.obj
            self.results['objvalue'] = bmwt
            table = {}
            n = 0  # tst
            if self.grp == 'gene':
                for gid in self.lst:
                    if self.mod.genes[gid].mnet==True:
                        n += 1
                        bds = self.pmod.getVarBounds(gid)
                        self.pmod.setVarBounds(gid, 0, self.eps) #set nonzero ub to avoid ill-conditioned
                        self.optsolveComb(solver=solver, *args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        table[gid] = {'st': st}
                        if st == 'optimal':
                            bm = self.pmod.getObjValue()
                            table[gid].update({'bm': str(bm)})
                            table[gid].update({'rate': str(float(bm)/bmwt)})
                            print str(n), gid, 'Growth='+table[gid]['bm']  # tst
                        else: print str(n), gid, 'Optimal_Not_Found'  # tst
                        self.pmod.setVarBounds(gid, bds[0], bds[1])

            elif self.grp == 'reaction':
                for rid in self.lst:
                    n += 1
                    bds = self.pmod.getVarBounds(rid)
                    self.pmod.setVarBounds(rid, 0, self.eps)
                    self.optsolveComb(solver=solver, *args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    table[rid] = {'st': st}
                    if st == 'optimal':
                        bm = self.pmod.getObjValue()
                        table[rid].update({'bm': str(bm)})
                        table[rid].update({'rate': str(float(bm)/bmwt)})
                        print str(n), rid, 'Growth=' + table[rid]['bm']  # tst
                    else: print str(n), rid, 'Optimal_Not_Found'  # tst
                    self.pmod.setVarBounds(rid, bds[0], bds[1])
            time_end = timer()
            self.results['runtime'] = time_end - time_start
            self.results['table'] = table
        else:
            self.results['objvalue'] = 'Optimal_Not_Found'
        #     raise myError, 'Cannot find optimal solution: '+ssbm+', '+ stbm


    def solveMedia(self, solver='glpk2', obj='', *args, **kwargs):
        time_start = timer()
        if self.prob and self.prob.media:
            outmd = {}  # multi-medium results
            for md in self.prob.media:
                print '>Medium: ' + md
                self.renewVarsBounds(self.prob.media[md])
                self.solve(solver=solver, obj=obj, *args, **kwargs)
                outmd[md] = cp.deepcopy(self.results)
            self.resultsMedia['outmd'] = outmd
            time_end = timer()
            self.resultsMedia['runtime'] = time_end - time_start
        else: self.solve(solver=solver, obj=obj, *args, **kwargs)



    #tst: for media-gene tests
    def solve2(self, biomass, *args, **kwargs):
        time_start = timer()
        self.setSolverOptions('glpk', 'tmlim', value=10)  # should be different option strategy for different solver
        self.setSolverOptions('glpk', 'dual', value='')
        self.setSolverOptions('glpk', 'nopresol', value='')
        self.setObjective(biomass)
        self.optsolve(*args, **kwargs)
        ssbm = str(self.pout.solver.status)
        stbm = str(self.pout.solver.termination_condition)
        bmwt = 1e-5
        if stbm == 'optimal':
            bmwt = self.pmod.getObjValue()
        self.results['solver_termination'] = stbm
        self.results['obj'] = self.obj
        self.results['objvalue'] = bmwt
        table = {}
        if self.grp == 'gene':
            X = 0 #tst
            for gid in self.lst:
                X += 1 #tst
                self.setSolverOptions('glpk', 'dual', value='')
                self.setSolverOptions('glpk', 'nopresol', value='')
                self.pmod.setVarBounds(gid, 0, 1e-8) #set ub=0 could be ill-conditioned
                self.optsolve(*args, **kwargs)
                st = str(self.pout.solver.termination_condition)
                if st != 'optimal':
                    self.options['glpk'].pop('nopresol')
                    self.optsolve(*args, **kwargs)
                    st = str(self.pout.solver.termination_condition)
                    if st != 'optimal':
                        self.setSolverOptions('glpk', 'nopresol', value='')
                        self.pmod.setVarBounds(gid, 0, 1e-7)
                        self.optsolve(*args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        if st != 'optimal':
                            self.options['glpk'].pop('dual')
                            self.pmod.setVarBounds(gid, 0, 1e-8)
                            self.optsolve(*args, **kwargs)
                            st = str(self.pout.solver.termination_condition)
                            if st != 'optimal':
                                self.pmod.setVarBounds(gid, 0, 1e-7)
                                self.optsolve(*args, **kwargs)
                                st = str(self.pout.solver.termination_condition)
                table[gid] = {'st': st}
                if st == 'optimal':
                    bm = self.pmod.getObjValue()
                    table[gid].update({'bm': str(bm)})
                    table[gid].update({'rate': str(float(bm)/bmwt)})
                self.pmod.setVarBounds(gid, 0, 1)
                print str(X), gid, table[gid] #tst
        time_end = timer()
        self.results['runtime'] = time_end - time_start
        self.results['table'] = table



    #tst: for testing media
    def solveProb(self, biomass='', prob = {}, *args, **kwargs):
        if (not biomass) and prob.obj: biomass = prob.obj
        if self.grp == 'gene' and prob.genes: self.lst = prob.genes
        elif self.grp == 'reaction' and prob.rxns: self.lst = prob.rxns
        if prob.constraints:
            for rid, bds in prob.constraints.items():
                self.pmod.setVarBounds(rid, bds[0], bds[1])
        if prob.media:
            cmed = {}
            for md in prob.media:
                print '-'*10+md+'-'*10 #tst
                for rid, bds in cmed.items(): #remove last medium
                    self.pmod.setVarBounds(rid, bds[0], bds[1])
                cmed = {}
                for rbd in prob.media[md]: #set new medium
                    rid = rbd[0]
                    cmed[rid] = self.pmod.vars[rid].bounds
                    self.pmod.setVarBounds(rid, rbd[1], rbd[2])
                self.solve2(biomass, *args, **kwargs)
                self.resultsMedia[md] = cp.deepcopy(self.results)
        else: self.solve2(biomass, *args, **kwargs)


    def write(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: KO,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.results['obj']]) + '\n')
        f.write('#TimeSpent: ' + str(self.results['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        table = self.results['table']
        if self.grp == 'gene':
            f.write('#Columns: ' + ', '.join(['ID', 'Growth', 'Ratio', 'Status']) + '\n')
            for id in table.keys():
                if table[id]['st'] == 'optimal':
                    f.write('\t'.join([id, table[id]['bm'], table[id]['rate'], table[id]['st']]) + '\n')
                else:
                    f.write('\t'.join([id, '', '', table[id]['st']]) + '\n')
        elif self.grp == 'reaction':
            f.write('#Columns: ' + ', '.join(['ID', 'Growth', 'Ratio', 'Status', 'Comment']) + '\n')
            for id in table.keys():
                if table[id]['st'] == 'optimal':
                    f.write('\t'.join([id, table[id]['bm'], table[id]['rate'], table[id]['st'], self.mod.getRxnComment(id)]) + '\n')
                else:
                    f.write('\t'.join([id, '', '', table[id]['st'], self.mod.getRxnComment(id)]) + '\n')
        f.close()


    def writeMedia(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: KO,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            f.write('-'*100 + '\n')
            f.write('#Medium: ' + md + '\n')
            f.write('#Objective: ' + ', '.join([str(results['objvalue']), results['solver_termination'], results['obj']]) + '\n')
            f.write('#TimeSpent: ' + str(results['runtime']) + ' (s)' + '\n')
            table = results['table']
            if self.grp == 'gene':
                f.write('#Columns: ' + ', '.join(['ID', 'Growth', 'Ratio', 'Status']) + '\n')
                for id in table.keys():
                    if table[id]['st'] == 'optimal':
                        f.write('\t'.join([id, table[id]['bm'], table[id]['rate'], table[id]['st']]) + '\n')
                    else:
                        f.write('\t'.join([id, '', '', table[id]['st']]) + '\n')
            elif self.grp == 'reaction':
                f.write('#Columns: ' + ', '.join(['ID', 'Growth', 'Ratio', 'Status', 'Comment']) + '\n')
                for id in table.keys():
                    if table[id]['st'] == 'optimal':
                        f.write('\t'.join([id, table[id]['bm'], table[id]['rate'], table[id]['st'], self.mod.getRxnComment(id)]) + '\n')
                    else:
                        f.write('\t'.join([id, '', '', table[id]['st'], self.mod.getRxnComment(id)]) + '\n')
        f.close()



class Tfba(Analysis):
    def __init__(self, mod, grp='reaction', lst='', *args, **kwargs):
        Analysis.__init__(self, mod, *args, **kwargs)
        self.grp = grp
        self.lst = lst
        if lst.__class__ == str: self.lst = lst.split()
        self.results = {}
        self.resultsMedia = {}  # results for each medium
        if self.prob: self.setProblem()
        if not self.lst:
            if self.grp == 'reaction': self.lst = self.pmod.rids
            elif self.grp == 'gene': self.lst = self.pmod.gids


    def setProblem(self):
        if self.prob.obj: self.setObjective(self.prob.obj)
        if self.prob.constraints: self.setVarsBounds(self.prob.constraints)
        if not self.lst:
            if self.grp == 'reaction':
                if self.prob.rxns:
                    self.lst = self.prob.rxns
                else: self.lst = self.pmod.rids
            elif self.grp == 'gene':
                if self.prob.genes:
                    self.lst = self.prob.genes
                else: self.lst = self.pmod.gids

    def solve(self, expr, array, solver='glpk2', med='', *args, **kwargs):
        time_start = timer()
        obj = {}
        for gid in self.pmod.gids:
            level = expr.level.get((array, gid), 0)
            if level!=0: obj.update({gid: level})

        self.setObjective(obj)
        self.optsolveComb(solver=solver, *args, **kwargs)
        self.results['solver_status'] = str(self.pout.solver.status)
        self.results['solver_termination'] = str(self.pout.solver.termination_condition)
        self.results['obj'] = self.obj

        if self.results['solver_termination'] == 'optimal':
            self.results['objvalue'] = self.pmod.getObjValue()
            self.results['table'] = self.pmod.getVarValues()
            if med:
                print med + '  ' + str(self.results['objvalue'])+'  '+self.results['solver_termination']
            else: print str(self.results['objvalue'])+'  '+self.results['solver_termination']
        else:
            if med:
                print med + '  ' + 'Optimal_Not_Found' + '  ' + self.results['solver_termination']
            else: print 'Optimal_Not_Found' + '  ' + self.results['solver_termination']
        time_end = timer()
        self.results['runtime'] = time_end - time_start


    def solveMedia(self, expr, array, solver='glpk2', *args, **kwargs):
        time_start = timer()
        if self.prob and self.prob.media:
            outmd = {}  # multi-medium results
            for md in self.prob.media:
                print '>Medium: '+md
                self.renewVarsBounds(self.prob.media[md])
                self.solve(expr, array, solver=solver, *args, **kwargs)
                outmd[md] = cp.deepcopy(self.results)
            self.resultsMedia['outmd'] = outmd
            time_end = timer()
            self.resultsMedia['runtime'] = time_end - time_start
        else: self.solve(expr, array, solver=solver, *args, **kwargs)

    def write_fba(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: TFBA,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.obj])+'\n')
        f.write('#TimeSpent: ' + str(self.results['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o+'='+str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'Value', 'Comment']) + '\n')
        table = self.results['table']
        if self.lst:
            for id in self.lst:
                f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(table[id]), self.mod.getRxnComment(id)])+'\n')
        else:
            for id, v in table.iteritems():
                f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(v), self.mod.getRxnComment(id)])+'\n')
        f.close()


    def writeMedia_fba(self, resultsfile):
        f = open(resultsfile, 'w')
        f.write('#Analysis: FBA,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            f.write('-' * 100 + '\n')
            f.write('#Medium: ' + md + '\n')
            f.write('#Objective: ' + ', '.join([str(results['objvalue']), results['solver_termination'], results['obj']]) + '\n')
            f.write('#TimeSpent: ' + str(results['runtime']) + ' (s)' + '\n')
            f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'Value', 'Comment']) + '\n')
            if results['solver_termination'] == 'optimal':
                table = results['table']
                if self.lst:
                    for id in self.lst:
                        f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(table[id]), self.mod.getRxnComment(id)]) + '\n')
                else:
                    for id, v in table.iteritems():
                        f.write('\t'.join([id, str(self.pmod.getVarLB(id)), str(self.pmod.getVarUB(id)), str(v), self.mod.getRxnComment(id)]) + '\n')
        f.close()



class Tfva(Analysis):
    #grp=reaction/gene
    def __init__(self, mod, grp='gene', lst='', *args, **kwargs):
        Analysis.__init__(self, mod, *args, **kwargs)
        self.grp = grp
        self.lst = lst
        if lst.__class__ == str: self.lst = lst.split()
        self.results = {}
        self.resultsMedia = {}  # results for each medium
        if self.prob: self.setProblem()
        if not self.lst:
            if self.grp == 'reaction': self.lst = self.pmod.rids
            elif self.grp == 'gene': self.lst = self.pmod.gids


    def setProblem(self):
        if self.prob.obj: self.setObjective(self.prob.obj)
        if self.prob.constraints: self.setVarsBounds(self.prob.constraints)
        if not self.lst:
            if self.grp == 'reaction':
                if self.prob.rxns:
                    self.lst = self.prob.rxns
                else: self.lst = self.pmod.rids
            elif self.grp == 'gene':
                if self.prob.genes:
                    self.lst = self.prob.genes
                else: self.lst = self.pmod.gids


    def solve(self, expr, array, solver='glpk2', *args, **kwargs):
        time_start = timer()
        obj = {}
        for gid in self.pmod.gids:
            level = expr.level.get((array, gid), 0)
            if level!=0: obj.update({gid: level})
        self.setObjective(obj)
        self.optsolveComb(solver=solver, *args, **kwargs)
        self.results['solver_status'] = str(self.pout.solver.status)
        self.results['solver_termination'] = str(self.pout.solver.termination_condition)
        self.results['obj'] = self.obj
        table = {}
        if self.results['solver_termination'] == 'optimal':
            self.results['objvalue'] = self.pmod.getObjValue()
            self.pmod.addConstraintObj(self.eps)
            if self.grp == 'gene':
                n = 0
                for id in self.lst:
                    if self.pmod.trvs[id].active == True:
                        self.setObjective(id, sense=1)
                        self.optsolveComb(solver=solver, *args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        table[id] = {'stmin': st}
                        table[id].update({'min': self.pmod.getVarValue(id)})
                        self.setObjective(id, sense=-1)
                        self.optsolveComb(solver=solver, *args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        table[id].update({'stmax': st})
                        table[id].update({'max': self.pmod.getVarValue(id)})
                        table[id].update({'state': str(expr.level.get((array, id), 0))})
                        n += 1
                        print str(n), id, 'min=' + str(table[id]['min']) + ':' + table[id]['stmin'], 'max=' + str(table[id]['max']) + ':' + table[id]['stmax']  # tst

            elif self.grp == 'reaction':
                n = 0
                for id in self.lst:
                    if self.pmod.vars[id].active == True:
                        self.setObjective(id, sense=1)
                        self.optsolveComb(solver=solver, *args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        table[id] = {'stmin': st}
                        table[id].update({'min': self.pmod.getVarValue(id)})
                        self.setObjective(id, sense=-1)
                        self.optsolveComb(solver=solver, *args, **kwargs)
                        st = str(self.pout.solver.termination_condition)
                        table[id].update({'stmax': st})
                        table[id].update({'max': self.pmod.getVarValue(id)})
                        n += 1
                        print str(n), id, 'min=' + str(table[id]['min']) + ':' + table[id]['stmin'], 'max=' + str(table[id]['max']) + ':' + table[id]['stmax']  # tst
            time_end = timer()
            self.results['runtime'] = time_end - time_start
            self.results['table'] = table
        else:
            self.results['objvalue'] = 'Optimal_Not_Found'
            print 'Optimal_Not_Found' + '  ' + self.results['solver_termination']


    def solveMedia(self, expr, array, solver='glpk2', *args, **kwargs):
        time_start = timer()
        if self.prob and self.prob.media:
            outmd = {}  # multi-medium results
            for md in self.prob.media:
                print '>Medium: '+md
                self.renewVarsBounds(self.prob.media[md])
                self.solve(expr, array, solver=solver, *args, **kwargs)
                outmd[md] = cp.deepcopy(self.results)
            self.resultsMedia['outmd'] = outmd
            time_end = timer()
            self.resultsMedia['runtime'] = time_end - time_start
        else: self.solve(expr, array, solver=solver, *args, **kwargs)


    def write(self, resultsfile):
        zcut = 1e-8 #zero cutoff for gene flux
        f = open(resultsfile, 'w')
        f.write('#Analysis: TFVA,  ' + self.cmt + '\n')
        f.write('#Objective: ' + ', '.join([str(self.results['objvalue']), self.results['solver_termination'], self.results['obj']]) + '\n')
        f.write('#TimeSpent: ' + str(self.results['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        table = self.results['table']
        if self.grp == 'gene':
            f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'State']) + '\n')
            for id in table.keys():
                if self.pmod.trvs[id].active == True:
                    activity = '0'
                    stmin = table[id]['stmin']
                    stmax = table[id]['stmax']
                    min = table[id]['min']
                    max = table[id]['max']
                    if stmin == 'optimal' and stmax == 'optimal':
                        if min > zcut:
                            activity = '1'
                        elif min <= zcut and max <= zcut:
                            activity = '-1'
                    f.write('\t'.join([id, str(self.pmod.trvs[id].lb), str(self.pmod.trvs[id].ub), str(min), stmin, str(max), stmax, activity, table[id]['state']]) + '\n')
        elif self.grp == 'reaction':
            f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'Comment']) + '\n')
            for id in table.keys():
                if self.pmod.vars[id].active == True:
                    activity = '0'
                    stmin = table[id]['stmin']
                    stmax = table[id]['stmax']
                    min = table[id]['min']
                    max = table[id]['max']
                    if stmin == 'optimal' and stmax == 'optimal':
                        if (min>zcut and max>zcut) or (min<-zcut and max<-zcut):
                            activity = '1'
                        elif abs(min) < zcut and abs(max) <= zcut:
                            activity = '-1'
                    f.write('\t'.join([id, str(self.pmod.vars[id].lb), str(self.pmod.vars[id].ub), str(min), stmin, str(max), stmax, activity, self.mod.getRxnComment(id)]) + '\n')
        f.close()


    def writeMedia(self, resultsfile):
        zcut = 1e-8 #zero cutoff for gene flux
        f = open(resultsfile, 'w')
        f.write('#Analysis: TFVA,  ' + self.cmt + '\n')
        f.write('#TimeSpent: ' + str(self.resultsMedia['runtime']) + ' (s)' + '\n')
        solopt = 'default'
        if self.options.get(self.solver, ''): solopt = ', '.join([o + '=' + str(v) for o, v in self.options[self.solver].items()])
        f.write('#Solver: ' + self.solver + ',  ' + 'Options: ' + solopt + '\n')
        for md, results in self.resultsMedia['outmd'].iteritems():
            f.write('-' * 100 + '\n')
            f.write('#Medium: ' + md + '\n')
            f.write('#Objective: ' + ', '.join([str(results['objvalue']), results['solver_termination'], results['obj']]) + '\n')
            f.write('#TimeSpent: ' + str(results['runtime']) + ' (s)' + '\n')
            table = results['table']
            if self.grp == 'gene':
                f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'State']) + '\n')
                for id in table.keys():
                    if self.pmod.trvs[id].active == True:
                        activity = '0'
                        stmin = table[id]['stmin']
                        stmax = table[id]['stmax']
                        min = table[id]['min']
                        max = table[id]['max']
                        if stmin == 'optimal' and stmax == 'optimal':
                            if min > zcut:
                                activity = '1'
                            elif min <= zcut and max <= zcut:
                                activity = '-1'
                        f.write('\t'.join([id, str(self.pmod.trvs[id].lb), str(self.pmod.trvs[id].ub), str(min), stmin, str(max), stmax, activity, table[id]['state']]) + '\n')
            elif self.grp == 'reaction':
                f.write('#Columns: ' + ', '.join(['ID', 'LB', 'UB', 'MIN', 'MIN_status', 'MAX', 'MAX_status', 'Activity', 'Comment']) + '\n')
                for id in table.keys():
                    if self.pmod.vars[id].active == True:
                        activity = '0'
                        stmin = table[id]['stmin']
                        stmax = table[id]['stmax']
                        min = table[id]['min']
                        max = table[id]['max']
                        if stmin == 'optimal' and stmax == 'optimal':
                            if (min>zcut and max>zcut) or (min<-zcut and max<-zcut):
                                activity = '1'
                            elif abs(min) < zcut and abs(max) <= zcut:
                                activity = '-1'
                        f.write('\t'.join([id, str(self.pmod.vars[id].lb), str(self.pmod.vars[id].ub), str(min), stmin, str(max), stmax, activity, self.mod.getRxnComment(id)]) + '\n')
        f.close()

