# TRN-FBA
A tool to perform constraint-based modelling on GSMN and TRMN models.

TRN-FBA is a constraint-based modelling (CBM) tool able to perform simulation for Transcriptional Regulatory Metabolic Network (TRMN) models, a proposed model that combines transcriptional regulatory network (TRN) and metabolic network, as well as Genome-scale Metabolic Network (GSMN) models. Six types of CBM functions were implemented including flux balance analysis, flux variability analysis, knock-out analysis, transcriptional flux balance analysis, and transcriptional flux variability analysis; both reactions and genes can be interrogated for all types of analyses. TRN-FBA was developed in Python 2.7 as a command-line tool, while you can also run it in Python console using its functions after importing it as a package. This tool can be used in any operating system such as Windows, Linux, and Mac system.

REQUIREMENTS

1)	Python 2.7+; we suggest installing Anaconda distribution of Python, which allows you easily install the following packages using command conda.

2)	Pyomo; it is a Python-based, open-source optimization modelling language with a diverse set of optimization capabilities. We applied Pyomo for building linear programming models and access the LP solvers, so far only GLPK and Gurobi solvers have been tested. Here is Pyomo link: https://pyomo.readthedocs.io/en/latest/index.html# .

3)	GLPK; GLPK solver, which can be installed by following this link: https://pyomo.readthedocs.io/en/latest/installation.html .

4)	Gurobi; alternatively you can use Gurobi solver which is faster than GLPK solver. To install Gurobi please go to this link: https://www.gurobi.com/documentation/8.1/quickstart_windows/installing_the_anaconda_py.html.


INSTALLATION

You just need to download package source codes to your specified directory, and add this directory to your system PATH, and then you can use it straight away.


COPYRIGHT
Copyright © 2020 Huihai Wu. This program is free software; you may redistribute it under the terms of the GNU General Public License v3.0. This program has no warranty.
