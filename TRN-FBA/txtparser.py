import component
from utility import *
import re



def parserEquationSide(expr, mod, ext):
    idcoef = {}
    mcs = expr.split(' + ')
    for mc in mcs:
        tt = mc.split()
        id = tt[0]
        coef =1
        if len(tt) == 2:
            id = tt[1]
            coef = float(tt[0])
        if not re.search(ext+'$', id):
            molecule = component.Molecule(id = id)
            mod.addMolecule(id, molecule)
            idcoef[id]=coef
    return idcoef


def parserEquation(expr, mod, ext):
    exprL, exprR = expr.split(' = ')
    lhs, rhs = '', ''
    lhs = parserEquationSide(exprL, mod, ext)
    if exprR:
        rhs = parserEquationSide(exprR, mod, ext)
    return lhs, rhs


def parserEruleGenes(rule, mod, rid):
    gids = re.findall(r'\b\w+\b', re.sub(r'AND|OR', '', rule))
    for gid in gids:
        if gid in mod.genes: mod.genes[gid].addRxn(rid)
        else:
            gene = component.Gene(id=gid)
            gene.addRxn(rid)
            gene.mnet = True
            mod.addGene(gid, gene)


def parserErule(rule, mod, ersn, nline):
    ruleExpr = rule
    id = ''
    ers = re.findall(r'(\([\w\s]*\))', rule)
    while ers:
        for er in ers:
            name = er.strip('(|)').strip()
            gates = re.findall(r' OR | AND ', name)
            if ('OR' in gates) and ('AND' in gates):
                raise ParserError(nline, ruleExpr)
            gate = gates[0].strip()
            sups = name.split(gates[0])
            sups = [x.strip() for x in sups]
            sname = gate+':'+''.join(sups)
            if sname in ersn: id = ersn[sname]
            else:
                id = 'E' + str(len(ersn))
                ersn[sname] = id
            erule = component.Erule(id=id, name=name, sname=sname, gate=gate, sups=sups)
            mod.addErule(id, erule)
            rule = rule.replace(er, id)
        ers = re.findall(r'(\([\w\s]*\))', rule)
    if re.search(' OR | AND ', rule):
        gates = re.findall(r' OR | AND ', rule)
        if ('OR' in gates) and ('AND' in gates):
            raise ParserError(nline, ruleExpr)
        gate = gates[0].strip()
        sups = rule.split(gates[0])
        sups = [x.strip() for x in sups]
        sname = gate + ':' + ''.join(sups)
        if sname in ersn:
            id = ersn[sname]
        else:
            id = 'E' + str(len(ersn))
            ersn[sname] = id
        erule = component.Erule(id=id, name=rule, sname=sname, gate=gate, sups=sups)
        mod.addErule(id, erule)
        return id
    else: return rule


def parserTranscription(expr):
    exprL, exprR = expr.split(' -> ')
    tfs = exprL.split(' + ')
    tgs = {} #target_gene:coef
    pgs = exprR.split(' + ')
    for pg in pgs:
        tt = pg.strip().split()
        id = tt[0]
        coef = 1
        if len(tt) == 2:
            id = tt[1]
            coef = float(tt[0])
        tgs[id] = coef
    return tfs, tgs


def parserTFconformation(expr):
    exprL, exprR = expr.split(' -> ')
    tt = exprL.split(' + ')
    tf = tt[0].strip()
    efs = tt[1]
    thr = 1e-3 #default thr
    if re.match(r'\d', efs):
        tt = efs.split(' ', 1)
        thr = float(tt[0])
        efs = tt[1]
    efs = efs.strip('(|)').strip().split(' OR ')
    tt = exprR.strip().split()
    tfe = tt[0]
    tfec = 1
    if len(tt)>1:
        tfe = tt[1]
        tfec = float(tt[0])
    return tf, efs, thr, tfe, tfec


def parserModelFile(modelfile, mod, ext):
    nline = 0
    ersn = {}  # erule_sname:eid
    for line in modelfile:
        nline += 1
        if re.match('#', line): continue
        tt = line.strip().split('\t')
        id = tt[0] #rxn id
        expr = tt[1]
        lb = float(tt[2])
        ub = float(tt[3])
        rule = tt[4].strip()
        cmt = tt[5]
        if re.match('^R', id):
            if not re.search(' = ', expr):
                raise ParserError(nline, line.strip())
            reaction = component.Reaction(id=id, equation=expr, lb=lb, ub=ub, rule=rule, cmt=cmt)
            lhs, rhs = parserEquation(expr, mod, ext)
            reaction.addLrhs(lhs, rhs)
            mod.addReaction(id, reaction)
            mod.updateCoef(id, lhs, 'l')
            mod.updateCoefm(id, lhs, 'l')
            if rhs:
                mod.updateCoef(id, rhs, 'r')
                mod.updateCoefm(id, rhs, 'r')
            if rule:
                parserEruleGenes(rule, mod, id)
                eid = parserErule(rule, mod, ersn, nline)
                mod.updateRxnErule(id, eid)

        elif re.match('^T', id) and rule in ['+', '-']:
            tfs, tgs = parserTranscription(expr)
            transcription = component.Transcription(id=id, expr=expr, tfs=tfs, tgs=tgs, type=rule)
            mod.addTranscription(id, transcription)
            mod.updateCoefg(tgs, id, rule)
            tfs = [x.rsplit('(', 1)[0] for x in tfs] #get gid for those in tfe format
            for gid in tfs: #add tf genes and target genes into mod.genes
                if gid in mod.genes:
                    mod.genes[gid].TF = True
                    mod.genes[gid].addTgs(tgs.keys())
                else:
                    gene = component.Gene(id=gid, TF=True)
                    gene.addTgs(tgs.keys())
                    mod.addGene(gid, gene)
            for gid in tgs.keys():
                if gid in mod.genes:
                    mod.genes[gid].addTfs(tfs)
                else:
                    gene = component.Gene(id=gid)
                    gene.addTfs(tfs)
                    mod.addGene(gid, gene)

        elif re.match('^T', id) and rule in ['e+', 'e-']:
            tf, efs, thr, tfe, tfec = parserTFconformation(expr)
            tfconformation = component.TFconformation(id=id, expr=expr, tf=tf, thr=thr, efs=efs, tfe=tfe, tfec=tfec, type=rule)
            mod.addTfconformations(id, tfconformation)
            if tf not in mod.genes:
                gene = component.Gene(id=tf, TF=True)
                mod.addGene(tf, gene)
        else:
            pass



def parserLinearExpr(obj):
    idcoef = {}
    mcs = obj.split(' + ')
    for mc in mcs:
        tt = mc.split()
        id = tt[0]
        coef = 1
        if len(tt) == 2:
            id = tt[1]
            coef = float(tt[0])
        idcoef[id] = coef
    return idcoef


def parserProbFile(probfile, prob):
    #nline = 0
    line = probfile.readline()
    while line:
        #nline += 1
        if re.match('#', line):
            line = probfile.readline()
            continue
        if re.match('>ANALYSIS', line):
            para = line.strip().split('\t')
            if len(para) > 1: prob.analysis = para[1].strip()
            line = probfile.readline()
        elif re.match('>OBJSENSE', line):
            para = line.strip().split('\t')
            if len(para) > 1: prob.sense = para[1].strip()
            line = probfile.readline()
        elif re.match('>OBJECTIVE', line):
            para = line.strip().split('\t')
            prob.obj = para[1].strip()
            #if len(para) > 1: prob.obj = parserLinearExpr(para[1].strip())
            line = probfile.readline()
        elif re.match('>REACTIONS', line):
            para = line.strip().split('\t')
            if len(para) > 1: prob.rxns = para[1].strip().split()
            line = probfile.readline()
        elif re.match('>GENES', line):
            para = line.strip().split('\t')
            if len(para) > 1: prob.genes = para[1].strip().split()
            line = probfile.readline()
        elif re.match('>CONSTRAINTS', line):
            line = probfile.readline()
            while line.strip():
                if re.match('>', line): break
                if re.match('#', line):
                    line = probfile.readline()
                    continue
                tt = line.strip().split('\t')
                prob.constraints[tt[0]] = (float(tt[1]), float(tt[2]))
                line = probfile.readline()
        elif re.match('>MEDIA', line):
            line = probfile.readline()
            while line.strip():
                if re.match('>', line): break
                if re.match('#', line):
                    line = probfile.readline()
                    continue
                tt = line.strip().split('\t')
                med = tt[0]
                if med in prob.media: prob.media[med].update({tt[1]: (float(tt[2]), float(tt[3]))})
                else: prob.media[med] = {tt[1]: (float(tt[2]), float(tt[3]))}
                line = probfile.readline()
        elif re.match('>SOLVER', line):
            para = line.strip().split('\t')
            if len(para) > 1: prob.solver = para[2].strip()
            line = probfile.readline()
            while line.strip():
                if re.match('>', line): break
                if re.match('#', line):
                    line = probfile.readline()
                    continue
                tt = line.strip().split('\t')
                prob.params[tt[0]] = tt[1]
                line = probfile.readline()
        else: line = probfile.readline()






def parserCompFile(componentfile, model):
    pass



def parserExprFile(exprfile, expr):
    line = exprfile.readline()
    tt = line.strip().split('\t')
    arrays = [x.strip() for x in tt[1:]]
    expr.arrays = arrays
    for line in exprfile:
        tt = line.strip().split('\t')
        gn = tt[0].strip()
        vs = [float(x) for x in tt[1:]] #correct
        for array, v in zip(arrays, vs):
            expr.level[array, gn] = v #correct
        expr.genes += [gn]


def parserSBMLmodelFile(exprfile, model):
    pass

