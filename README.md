# CVRP-Formulations

The cvrp module contains the CVRP class that is designed to
handle and solve instances of the CVRP obtained from the CVRPLIB,
using the DOcplex interface of the CPLEX solver.
It depends on the docplex, matplotlib, seaborn, numpy, urllib, and
scipy packages.

It can handle instances with euclidean two-dimensional distance 
edge weight types.
The class can be instanced using an url as follows:

import cvrp
instance = cvrp.CVRP.from_url('http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n16-k8.vrp')

The class implements 4 formulations (MTZ-L, GG, BHM, MCF) that can be
accesed through the following methods:

-.mtz()

-.gg()

-.bhm()

-.mcf()


These methods returns a docplex.mp.model.Model instance.

Alternatively, the .solve() method can be used to solve the instances
by passing one of the formulation names ('mtz', 'gg', 'bhm','mcf') 
through the 'formulation' parameter. It gives an ouput file and a figure
representation of the solution found. Also, a timelimit parameter can be
passed.

Example:
instance.solve(formulation='bhm', timelimit=2000)

The 'solve_instances' contains an example of solving 10 instances of the 
CVRPLIB using the .solve() method. It requires an internet connection to
work.
