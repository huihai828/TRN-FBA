
class Molecule:
    def __init__(self, id='', name='', formula='', charge=''):
        self.id = id
        self.name = name
        self.formula = formula
        self.charge = charge



class Reaction:
    def __init__(self, id='', name='', subsys='', equation='', lb=0, ub=1000, rule='', cmt=''):
        self.id = id
        self.name = name
        self.subsys = subsys
        self.equation = equation
        self.lhs = {} #molecule_id : coef
        self.rhs = {}
        self.rule = rule
        self.ruleid = ''
        #self.rule = dict() #AND/OR tree, e.g. ( A AND B ) OR ( C AND D ), OR1:{AND1:(A,B), AND2:{C,D)}
        self.reversible = False
        self.lb = lb
        self.ub = ub
        self.EC = ''
        self.cmt = cmt


    def addLrhs(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs


    def addBounds(self, lb, ub):
        self.lb = lb
        self.ub = ub


    def getBounds(self):
        return (self.lb, self.ub)


    def getComment(self):
        return self.cmt


    def setRuleid(self, eid):
        self.ruleid = eid


    def show(self):
        print 'ID: ' + self.id
        print 'Equation: ' + self.equation
        print 'LB: ' + str(self.lb)
        print 'UB: ' + str(self.ub)
        print 'Rule: ' + self.rule
        print 'Comment: ' + self.cmt



class Gene:
    def __init__(self, id='', name='', locus='', product='', TF=False):
        self.id = id
        self.name = name
        self.locus = locus
        self.product = product
        self.rxns = [] #rxns ruled by this enzyme gene
        self.TF = TF
        self.tgs = [] #target genes if it is TF
        self.tfs = [] #TF genes on this gene
        self.mnet = False


    def addRxn(self, rid):
        self.rxns += [rid]


    def addTgs(self, tgs):
        for tg in tgs:
            if tg not in self.tgs:
                self.tgs += [tg]


    def addTfs(self, tfs):
        for tf in tfs:
            if tf not in self.tfs:
                self.tfs += [tf]


    def show(self):
        print 'ID: ' + self.id
        print 'Metnet: ' + str(self.mnet)
        print 'TF: ' + str(self.TF)
        print 'Reactions: ' + ' '.join(self.rxns)
        print 'TFgenes: ' + ' '.join(self.tfs)
        print 'TargetGenes: ' + ' '.join(self.tgs)



class Erule: #enzyme complex gene rule
    def __init__(self, id='', name='', sname='', gate='', sups=[]):
        self.id = id #erule id, e.g. E123
        self.name = name #erule expression, e.g. ( A AND B ) OR ( C AND D )
        self.sname = sname #short hand name, e.g. AND:g1g2g3
        self.gate = gate #AND, OR
        self.sups = sups #super erule ids


    def show(self):
        print 'ID: ' + self.id
        print 'Name: ' + self.name
        print 'Sname: ' + self.sname
        print 'Gate: ' + self.gate
        print 'Supers: ' + ' '.join(self.sups)



class Transcription:
    def __init__(self, id='', name='', expr='', tfs=[], tgs={}, type=''):
        self.id = id
        self.name = name
        self.expr = expr
        self.tfs = tfs #list of input TFs
        self.tgs = tgs #list of target genes with coefs of probability e.g. p(A|tf1 and tf2)
        self.type = type #+: productive, -:repressive, +-: dual,
        self.lb = 0
        self.ub = 1
        self.cmt = ''


    def show(self):
        print 'ID: ' + self.id
        print 'Expression: ' + self.expr
        print 'TFs: ' + ' '.join(self.tfs)
        print 'Targets: ' + ' '.join(self.tgs)
        print 'Type: ' + self.type



class TFconformation:
    def __init__(self, id='', name='', expr='', tf='', efs=[], thr=0.0001, tfe=[], tfec=1, type=''):
        self.id = id
        self.name = name
        self.tf = tf
        self.efs = efs #list of ecffectors (OR)
        self.thr = thr #threshold or efs to be fully active
        self.tfe = tfe #TF-effector
        self.tfec = tfec #tfe coef: probability of conformation
        self.expr = expr
        self.type = type #e+:activate TF, e-:inactivate TF
        self.lb = 0
        self.ub = 1
        self.cmt = ''


    def show(self):
        print 'ID: ' + self.id
        print 'Expression: ' + self.expr
        print 'TF: ' + self.tf
        print 'Effectors: ' + ' '.join(self.efs)
        print 'TF_effectors: ' + self.tfe +': '+str(self.tfec)
        print 'Threshold: ' + str(self.thr)
        print 'Type: ' + self.type




