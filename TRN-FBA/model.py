import component

class Model:
    def __init__(self, id='', name=''):
        self.id = id
        self.name = name
        self.molecules = {}
        self.coef ={} #(mid, rid):coef
        self.coefm ={} #mid:{rid:coef}
        self.genes = {}
        self.reactions = {} #rid:reaction
        self.transcriptions = {} #tra
        self.tfconformations = {} #tfc
        self.compartments = {} #cid:name
        self.erules = {}
        self.rxner = {} #rxn_id: its erule_id
        self.coefg = {} #(tgid, '+/-'):{trid:coef}; +:productive, -:respressive
        #self.bmr = '' #biomass rxn


    def addMolecule(self, id, molecule):
        self.molecules.update({id:molecule})


    def addReaction(self, id, reaction):
        self.reactions.update({id:reaction})


    def addErule(self, id, erule):
        self.erules.update({id:erule})


    def addGene(self, id, gene):
        self.genes.update({id:gene})


    def addTranscription(self, id, transcription):
        self.transcriptions.update({id:transcription})


    def addTfconformations(self, id, tfconformation):
        self.tfconformations.update({id:tfconformation})


    def updateCoef(self, rid, hs, lr ):
        for mid, c in hs.iteritems():
            self.coef[mid, rid] = c * {'l': -1, 'r': 1}[lr]


    def updateCoefm(self, rid, hs, lr ):
        for mid, c in hs.iteritems():
            if mid in self.coefm:
                self.coefm[mid].update({rid: c * {'l': -1, 'r': 1}[lr]})
            else:
                self.coefm[mid] = {}
                self.coefm[mid].update({rid: c * {'l': -1, 'r': 1}[lr]})


    def updateRxnErule(self, rid, eid):
        self.rxner.update({rid:eid})
        self.reactions[rid].setRuleid(eid)


    def updateCoefg(self, tgcoef, tid, type):
        for gid, coef in tgcoef.items():
            if (gid, type) in self.coefg:
                self.coefg[gid, type].update({tid:coef})
            else:
                self.coefg[gid, type] = {tid:coef}


    def mol2RxnCoefs(self, mid):
        return self.coefm[mid]


    def getRxnComment(self, id):
        if id in self.reactions:
            return self.reactions[id].getComment()
        else: return ''

    def isTrmn(self): #indicate if model is GSMN or TRMN model
        if self.transcriptions: return True
        else: return False




class Problem:
    def __init__(self):
        self.analysis = ''
        self.sense = ''
        self.obj = '' #obj expr
        self.rxns = []
        self.genes = []
        self.constraints = {} #rxn/expr:(lb, ub)
        self.media = {} #medium_id: [(rxn1, lb, ub), (rxn2, lb, ub)]
        self.solver = ''
        self.options = {} #option:value



class ExpressionData:
    def __init__(self):
        self.arrays = [] #list of array ids (col)
        self.genes = [] #list of gene ids (raw)
        self.level = {} #discretized expr levels (array, gene) -> expr
