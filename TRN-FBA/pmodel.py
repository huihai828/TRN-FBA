import pyomo.environ as pe
import model
import component
import re

class Pmodel(pe.ConcreteModel):
    def __init__(self, mod='', *args, **kwargs):
        pe.ConcreteModel.__init__(self, *args, **kwargs)
        #self.name = name
        self.buildPmodel(mod)

    def buildPmodel(self, mod):
        #build pmodel for reactions part
        self.rids = mod.reactions.keys() #rxn id list
        self.gids = mod.genes.keys()  # gene id list
        mids = mod.coefm.keys()
        def getVarBds(pmod, i):
            return mod.reactions[i].getBounds()
        self.rset = pe.Set(initialize = self.rids)
        self.vars = pe.Var(self.rset, bounds=getVarBds) #vars for rxn flux
        def getCstExpr(pmod, i):
            return pe.quicksum(c*self.vars[r] for r, c in mod.mol2RxnCoefs(i).iteritems()) == 0
        self.mset = pe.Set(initialize=mids)
        self.csts = pe.Constraint(self.mset, rule=getCstExpr) #csts for rxn related constraints

        # build pmodel for transcriptions part
        if mod.isTrmn(): #for TRMN model
            eids = mod.erules.keys() #erule id
            tids = mod.transcriptions.keys() #transcription id
            cids = mod.tfconformations.keys() #tfconfromation id
            tfes = [mod.tfconformations[cid].tfe for cid in cids] #TFeffector id
            self.tset = pe.Set(initialize=self.gids+eids+tids+cids+tfes)
            self.trvs = pe.Var(self.tset, bounds=(0,1))
            self.cler = pe.ConstraintList()
            for eid in eids: #set erule constraints
                gate = mod.erules[eid].gate
                sups = mod.erules[eid].sups
                if gate == 'AND':
                    for sup in sups:
                        self.cler.add(self.trvs[eid] <= self.trvs[sup])
                elif gate == 'OR':
                    self.cler.add(self.trvs[eid] <= pe.quicksum(self.trvs[i] for i in sups))
            for rid, eid in mod.rxner.iteritems(): #set reaction constraints with its erule
                Vmax = mod.reactions[rid].ub
                Vmin = mod.reactions[rid].lb
                if Vmax > 0:
                    self.cler.add(self.vars[rid] <= Vmax*self.trvs[eid])
                if Vmin < 0:
                    self.cler.add(self.vars[rid] >= Vmin*self.trvs[eid])
            self.clgt = pe.ConstraintList() #set constrainsts between genes and transcriptions
            for gtype, tc in mod.coefg.iteritems():
                gid = gtype[0]
                type = gtype[1]
                if type == '+': #productive constraints
                    if len(tc) == 1:
                        self.clgt.add(self.trvs[gid] <= pe.quicksum(c * self.trvs[t] for t, c in tc.items()))
                    else:
                        self.clgt.add(self.trvs[gid] <= 0.8 * pe.quicksum(c * self.trvs[t] for t, c in tc.items()))
                elif type == '-': #repressive constraints
                    for t, c in tc.items():
                        self.clgt.add(self.trvs[gid] <= 1 - c * self.trvs[t])
            for tid in tids: #set transcription constraints by its tfs
                tfs = mod.transcriptions[tid].tfs
                type = mod.transcriptions[tid].type
                if type == '+':
                    for tf in tfs: self.clgt.add(self.trvs[tid] <= self.trvs[tf])
                elif type == '-': #for repressive transcription, '<=' means no constraint, so should use '=='
                    for tf in tfs: self.clgt.add(self.trvs[tid] == self.trvs[tf])
            self.clef = pe.ConstraintList()  # set constrainsts between tfconformation and effector
            for cid in cids:
                tf = mod.tfconformations[cid].tf
                efs = mod.tfconformations[cid].efs
                thr = mod.tfconformations[cid].thr
                fthr = float(1)/thr
                tfe = mod.tfconformations[cid].tfe
                tfec = mod.tfconformations[cid].tfec
                type = mod.tfconformations[cid].type
                efsrxn = []
                for ef in efs:
                    Ref = 'R_EF_' + ef #demand rxn id
                    efsrxn += [Ref]
                    self.setVarBounds(Ref, 0, thr)
                self.clef.add(self.trvs[cid] <= self.trvs[tf])
                if type=='e+':
                    self.clef.add(self.trvs[cid] <= fthr*pe.quicksum(self.vars[r] for r in efsrxn))
                elif type=='e-':
                    for r in efsrxn:
                        self.clef.add(self.trvs[cid] <= 1-fthr*self.vars[r])
                self.clef.add(self.trvs[tfe] <= tfec*self.trvs[cid])

        else: #for GSMN model
            eids = mod.erules.keys()  # erule id
            self.tset = pe.Set(initialize=self.gids + eids)
            self.trvs = pe.Var(self.tset, bounds=(0, 1))
            self.cler = pe.ConstraintList()
            for eid in eids:  # set erule constraints
                gate = mod.erules[eid].gate
                sups = mod.erules[eid].sups
                if gate == 'AND':
                    for sup in sups:
                        self.cler.add(self.trvs[eid] <= self.trvs[sup])
                elif gate == 'OR':
                    self.cler.add(self.trvs[eid] <= pe.quicksum(self.trvs[i] for i in sups))
            for rid, eid in mod.rxner.iteritems():  # set reaction constraints with its erule
                Vmax = mod.reactions[rid].ub
                Vmin = mod.reactions[rid].lb
                if Vmax > 0:
                    self.cler.add(self.vars[rid] <= Vmax * self.trvs[eid])
                if Vmin < 0:
                    self.cler.add(self.vars[rid] >= Vmin * self.trvs[eid])




    def addObjective(self, obj, sense, mod):
        self.pobj = {}
        for id in obj:
            if (id in self.vars) or (id in self.trvs): self.pobj[id] = obj[id]
            elif id in self.mset: # !wrong: for metabolie, set a demand exchange rxn for this met
                pass  # add var of this demand rxn and change cst of this met
                # for rid, c in mod.mol2RxnCoefs(id).items():
                #     if rid in self.rset and c > 0: self.pobj[rid] = c
        if hasattr(self, 'obj'): self.del_component(self.obj)
        #self.objexpr = pe.quicksum(c * self.vars[r] for r, c in self.pobj.items())
        self.objexpr = 0
        for id, c in self.pobj.items():
            if id in self.vars: self.objexpr = self.objexpr + c * self.vars[id]
            elif id in self.trvs: self.objexpr = self.objexpr + c * self.trvs[id]
        self.obj = pe.Objective(expr=self.objexpr, sense=sense)


    def setVarBounds(self, id, lb=float('-Inf'), ub=float('Inf')):
        if id in self.vars:
            self.vars[id].setlb(lb)
            self.vars[id].setub(ub)
        elif id in self.trvs:
            self.trvs[id].setlb(lb)
            self.trvs[id].setub(ub)


    def addConstraint(self, id, expr, lb=float('-Inf'), ub=float('Inf')):
        #if id in self.csts: self.del_component(self.csts[id])
        if id not in self.csts: self.mset.add(id)
        #self.mset.add(id)
        self.csts.add(id, expr= pe.inequality(lb, pe.quicksum(c * self.vars[r] for r, c in expr.items()), ub))


    #eps: tolerance coef for fixing the vars
    def addConstraintObj(self, eps):
        #if '_OBJ_' in self.csts: self.del_component(self.csts['_OBJ_']) #will delete whole csts
        if '_OBJ_' not in self.csts: self.mset.add('_OBJ_')
        objvalue = self.getObjValue()
        lb = objvalue - eps
        ub = objvalue + eps
        self.csts.add('_OBJ_', expr=pe.inequality(lb, self.objexpr, ub))


    def solve(self, solver, options=None, *args, **kwargs):
        optimizer = pe.SolverFactory(solver)
        if options:
            for k, v in options.items():
                optimizer.options[k] = v
        results = optimizer.solve(self, *args, **kwargs)
        #results = optimizer.solve(self, options=options, *args, **kwargs)
        return results


    #set stale to avoid ValueError: No value for uninitialized NumericValue object trvs[sodC]
    def getVarValues(self):
        varval = {}
        for v in self.component_data_objects(pe.Var, active=True):
            if not v.stale:
                id = re.findall(r'\w+\[(.+)\]', pe.name(v))[0]
                varval[id] = pe.value(v)
        return varval


    def getVarValue(self, id):
        if id in self.vars:
            return pe.value(self.vars[id])
        else: return pe.value(self.trvs[id])


    #get variable bounds from pmodel, return (lb, ub)
    def getVarBounds(self, id):
        if id in self.vars:
            return self.vars[id].bounds
        elif id in self.trvs:
            return self.trvs[id].bounds


    #get variable lower bound
    def getVarLB(self, id):
        if id in self.vars:
            return self.vars[id].lb
        elif id in self.trvs:
            return self.trvs[id].lb


    # get variable upper bound
    def getVarUB(self, id):
        if id in self.vars:
            return self.vars[id].ub
        elif id in self.trvs:
            return self.trvs[id].ub


    def getObjValue(self):
        return pe.value(self.obj)


    # get variable upper bound
    def getVarActive(self, id):
        if id in self.vars:
            return self.vars[id].active
        elif id in self.trvs:
            return self.trvs[id].active


    #print individual component info.
    def myprint(self, id):
        if id in self.vars:
            print '\t'.join([str(self.vars[id]), str(self.vars[id].lb), str(self.vars[id].ub), str(self.vars[id].value)])
        elif id in self.trvs:
            print '\t'.join([str(self.trvs[id]), str(self.trvs[id].lb), str(self.trvs[id].ub), str(self.trvs[id].value)])
        elif id in self.csts:
            print '\t'.join([str(self.csts[id].body), str(self.csts[id].lower), str(self.csts[id].upper)])


    def writeModel(self, outfile, mod):
        with open(outfile, 'w') as f:
            f.write('CONSTRAINTS: Metabolism' + '\n')
            for id in self.mset:
                f.write('\t'.join([id, str(self.csts[id].body), str(self.csts[id].lower), str(self.csts[id].upper)])+'\n')
            f.write('\n' + 'VARIABLES: Metabolism' + '\n')
            for id in self.rset:
                f.write('\t'.join([str(self.vars[id]), str(self.vars[id].lb), str(self.vars[id].ub)])+'\n')
            if mod.isTrmn():
                f.write('\n' + '-'*80 + 'TR_part' + '-'*80 + '\n')
                f.write('CONSTRAINTS: Enzyme_rule' + '\n')
                for i in self.cler:
                    f.write('\t'.join([str(i), str(self.cler[i].body), str(self.cler[i].lower), str(self.cler[i].upper)]) + '\n')
                f.write('\n' + 'CONSTRAINTS: Transcription' + '\n')
                for i in self.clgt:
                    f.write('\t'.join([str(i), str(self.clgt[i].body), str(self.clgt[i].lower), str(self.clgt[i].upper)]) + '\n')
                f.write('\n' + 'CONSTRAINTS: TFconformation' + '\n')
                for i in self.clef:
                    f.write('\t'.join([str(i), str(self.clef[i].body), str(self.clef[i].lower), str(self.clef[i].upper)]) + '\n')
                f.write('\n' + 'VARIABLES: TR_part' + '\n')
                for id in self.tset:
                    f.write('\t'.join([str(self.trvs[id]), str(self.trvs[id].lb), str(self.trvs[id].ub)]) + '\n')





