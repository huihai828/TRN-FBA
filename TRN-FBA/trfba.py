import argparse
import utility
import analysis as ta
import re


"""
#command line
-m --input input sfba model file, SBML file
-c --compo sfba model components table file
-e --expr gene expression file
-o --obj objective expression, e.g max:R_biomass, "min:-0.2 R1 + 0.5 R2"
-a --analysis choose analysis (obj, fba, fva, ko)
-l --list choose list of reactions/genes for analysis
-s --solver choose slover (glpk, gurobi) can followed by slover params (e.g. "glpk exact tmlim:6 tol_bnd:1e-6")
-p --pfile problem file (will be overwriten if in command line)
-f --outfile output results file name
-x --extag boudary/exteranl metabolite tag
--wpmod write out pyomo model
--wsfba export sfba model file
--wcomp export sfba component file
--wsbml export SBML model file
"""



def argsParser():
    parser = argparse.ArgumentParser(description='TRMN FBA analysis')
    parser.add_argument('-m', '--model', required=True,
                        help='input model file in SurreyFBA format')
    #parser.add_argument('-mc', '--compo',
    #                    help='input model components table file')
    parser.add_argument('-e', '--expr',
                        help='input gene expression file')
    parser.add_argument('-a', '--analysis', required=True, choices=['obj', 'fba', 'fvaR', 'fvaG', 'koR', 'koG', 'tfbaR', 'tfbaG', 'tfvaR', 'tfvaG'],
                        help='analysis types: obj - find optimized objective, \
                             fba - find a flux balance solution, \
                             fvaG/fvaR - flux variability analysis for genes/reactions, \
                             koG/koR - knockout analysis for genes/reactions, \
                             tfvaG/tfvaR - transcriptional flux variability analysis for genes/reactions')
    parser.add_argument('-o', '--obj',
                        help='objective function or expression array ID')
    parser.add_argument('-l', '--lst',
                        help='list of reactions/genes for analysis')
    parser.add_argument('-s', '--solver', default='glpk2',
                        help='choose slover such as glpk, gurobi, etc. default: glpk2 (glpk with choosen options)')
    parser.add_argument('-i', '--opt',
                        help='options for the slover e.g. -i "dual nopresol tmlim:6"')
    parser.add_argument('-p', '--pfile',
                        help='problem file can be overwriten by those in command line')
    parser.add_argument('-f', '--outfile',
                        help='output results file name')
    parser.add_argument('-x', '--extag', default='_b',
                        help='boudary/exteranl metabolite suffix')
    # parser.add_argument('-r', '--array',
    #                     help='specify array ID of expression data for tfva')
    parser.add_argument('-c', '--csts', action='append',
                        help='add variable constraint (id lb ub) e.g. -c "R1 1.5 10"')
    #parser.add_argument('-g', '--msg',
    #                    help='standard output message level: 0 - no message; 1 - process bar; \
    #                    2 - list output; 3 - solver info and list output')
    parser.add_argument('-t', '--eps', type=float, default=0,
                        help='a constraints fault tolerance for fva and tfva, default: 0')

    args = parser.parse_args()
    return args



def parseSolverOption(argopt):
    optvals = {}
    tt = argopt.strip().split()
    for optval in tt:
        optval = optval.strip().split(':')
        opt = optval[0]
        val = ''
        if len(optval)==2: val = optval[1]
        optvals.update({opt: val})
    return optvals



def main():
    args = argsParser()
    model = ta.readModel(args.model, name=args.model)
    #ta.setTFeffectorID(model)
    argcst = {}
    if args.csts:
        for cst in args.csts:
            cst = cst.strip().split()
            argcst[cst[0]] = (float(cst[1]), float(cst[2]))
    prob = ''
    if args.pfile: prob = ta.readProblem(args.pfile)
    options = {}
    if args.opt: options[re.sub('\d$', '', args.solver)] = parseSolverOption(args.opt)
    grp = {'G': 'gene', 'R': 'reaction'}.get(args.analysis[-1], '')

    if args.analysis in ['obj', 'fba']:
        fba = ta.Fba(model, lst=args.lst, options=options, prob=prob, cmt='Model: '+args.model)
        if argcst: fba.setVarsBounds(argcst)
        fba.solveMedia(solver=args.solver, obj=args.obj)
        if args.analysis=='obj':
            if args.outfile:
                if prob and prob.media: fba.writeMedia_obj(args.outfile)
                else: fba.write_obj(args.outfile)
        elif args.analysis=='fba':
            if args.outfile:
                if prob and prob.media: fba.writeMedia_fba(args.outfile)
                else: fba.write_fba(args.outfile)

    elif re.match('tfba', args.analysis):
        if not args.expr: raise Exception, 'ArgError: need expression file'
        expr = ta.readExprData(args.expr)
        fba = ta.Tfba(model, grp=grp, lst=args.lst, eps=args.eps, options=options, prob=prob, cmt='Model: '+args.model)
        if argcst: fba.setVarsBounds(argcst)
        fba.solveMedia(expr, args.obj, solver=args.solver)
        if args.outfile:
            if prob and prob.media: fba.writeMedia_fba(args.outfile)
            else: fba.write_fba(args.outfile)

    elif re.match('fva', args.analysis):
        fva = ta.Fva(model, grp=grp, lst=args.lst, eps=args.eps, options=options, prob=prob, cmt='Model: '+args.model)
        if argcst: fva.setVarsBounds(argcst)
        fva.solveMedia(solver=args.solver, obj=args.obj)
        if args.outfile:
            if prob and prob.media:
                fva.writeMedia(args.outfile)
            else:
                fva.write(args.outfile)

    elif re.match('ko', args.analysis):
        ko = ta.Ko(model, grp=grp, lst=args.lst, eps=args.eps, options=options, prob=prob, cmt='Model: '+args.model)
        if argcst: ko.setVarsBounds(argcst)
        ko.solveMedia(solver=args.solver, obj=args.obj)
        if args.outfile:
            if prob and prob.media:
                ko.writeMedia(args.outfile)
            else:
                ko.write(args.outfile)

    elif re.match('tfva', args.analysis):
        if not args.expr: raise Exception, 'ArgError: need expression file'
        expr = ta.readExprData(args.expr)
        tfva = ta.Tfva(model, grp=grp, lst=args.lst, eps=args.eps, options=options, prob=prob, cmt='Model: '+args.model)
        if argcst: tfva.setVarsBounds(argcst)
        tfva.solveMedia(expr, args.obj, solver=args.solver)
        if args.outfile:
            if prob and prob.media:
                tfva.writeMedia(args.outfile)
            else:
                tfva.write(args.outfile)



if __name__ == '__main__':
    try:
        main()

    except utility.ParserError, e:
        print 'ParserError at line '+e.line+ ': '+e.msg
    # except Exception, e:
    #    print e
    # finally:
    #     pass #close files
        #print 'finally run'
