import cvrp

set_p = ['http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n16-k8.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n19-k2.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n20-k2.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n23-k8.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n40-k5.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/A/A-n44-k6.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/A/A-n46-k7.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/A/A-n32-k5.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/A/A-n37-k6.vrp',
            'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/A/A-n62-k8.vrp']

formulations = ['mtz', 'gg', 'bhm', 'mcf']

for url in set_p:
    instance = cvrp.CVRP.from_url(url)
    for formulation in formulations:
        sol = instance.solve(formulation=formulation, timelimit=2000)