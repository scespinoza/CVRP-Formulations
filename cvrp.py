import urllib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns; sns.set()
from docplex.cp.parameters import CpoParameters
from docplex.mp.model import Model
from docplex.mp.conflict_refiner import ConflictRefiner
from scipy.spatial import distance_matrix

def read_cvrp(url, to_dat=False):

    cvrp = urllib.request.urlopen(url)
    header = {}
    data = []

    for line in cvrp:
        thisline = line.decode('UTF-8').split(':')
        if len(thisline) > 1:
            header[thisline[0].strip()] = ':'.join(thisline[1:]).strip()
        else:
            data.append(line.decode('UTF-8').strip())
   

    n = int(header['DIMENSION'])
    Q = int(header['CAPACITY'])
    comment = header['COMMENT'].split(':')
    K = int(comment[1].strip()[0])
    optimal = int(comment[2].strip().replace(')', ''))
    
    start_nodes = data.index('NODE_COORD_SECTION')
    nodes = ' '.join(data[start_nodes + 1:start_nodes + n + 1])
    nodes = np.fromstring(nodes, sep=' ').reshape(n, 3)[:, 1:]

    start_demand = data.index('DEMAND_SECTION')
    demands = ' '.join(data[start_demand + 1: start_demand + n + 1])
    demands = np.fromstring(demands, sep=' ').reshape(n, 2)[:, 1]

    start_depot = data.index('DEPOT_SECTION')
    depot = ' '.join(data[start_depot + 1:start_depot + 3])
    depot = np.fromstring(depot, sep=' ').reshape(1, 2)
    
    distances = np.rint(distance_matrix(nodes, nodes))
    np.fill_diagonal(distances, 9999)

    if to_dat:
        file = open(url.split('/')[-1].split('.')[0] + '.dat', 'w')
        file.write('n =' + str(n) + ';\n')
        file.write('K =' + str(K) + ';\n')
        file.write('Q =' + str(Q) + ';\n')
        file.write('q = [')
        for i in range(n):
                file.write(str(demands[i]) + ' ')
        file.write('];\n')
        file.write('d = [')
        for i in range(n):
            for j in range(n):
                file.write(str(distances[i, j]) + ' ')
            file.write('\n')
        file.write('];')
        file.close()

    return header, n, K, Q, demands, distances, nodes, optimal




class CVRP():

    def __init__(self, n, K, Q, demands, distances, description=None, nodes=None, optimal=None):
        self.n = n
        self.K = K
        self.Q = Q
        self.q = demands
        self.nodes = nodes
        self.d = distances
        self.description = description
        self.current_sol = None
        self.current_formulation = None
        self.current_obj = None
        self.optimal = optimal

    def get_routes(self):

        def get_next(node):
            for j in range(0, self.n):
                if node != j:
                    if self.current_sol[node, j].solution_value == 1:
                        return j

        n_routes = 0
        routes = {}
        root = 0

        for j in range(1, self.n):
            if self.current_sol[root, j].solution_value == 1:
                routes[n_routes] = [0, j]
                n_routes += 1
        
        for i in range(n_routes):
            current_node = routes[i][-1]
            while current_node!=0:
                routes[i].append(get_next(current_node))
                current_node = routes[i][-1]

        return routes
    
    @classmethod
    def from_url(cls, url):
        
        description, n, K, Q, demands, distances, nodes, optimal = read_cvrp(url, to_dat=False)   
        return cls(n, K, Q, demands, distances, description, nodes, optimal)

    def __str__(self):
        string = ''
        for key in self.description:
            string += key + ': ' + self.description[key] + ('\n')
        return string
    
    def plot(self, sol=True):
        fig = plt.figure()
        if not sol:
            plt.scatter(self.nodes[1:, 0], self.nodes[1:, 1])
            plt.scatter(self.nodes[0, 0], self.nodes[0, 1])
            plt.show()
        else:
            routes = self.get_routes()
            n_routes = max([i for i in routes]) + 1
            cmap = cm.get_cmap('Set1', n_routes)
            
                        
            if self.current_formulation == 'BHM':
                for i in range(self.n):
                  if int(self.current_sol[i, self.n].solution_value) == 1 or int(self.current_sol[self.n, i].solution_value):
                      plt.plot(self.nodes[[0, i], 0], self.nodes[[0, i], 1], lw=2, color='red', zorder=1, alpha=0.8)
                for i in range(self.n):
                    for j in range(self.n):
                        if i!=j:
                            if int(self.current_sol[i, j].solution_value) == 1:
                                plt.plot(self.nodes[[i, j], 0], self.nodes[[i, j], 1], lw=2, color='red', zorder=1, alpha=0.8)
            else:
                for i in routes:
                    for j in range(len(routes[i])-1):
                        src, dest = routes[i][j:j+2]
                        plt.plot(self.nodes[[src, dest], 0], self.nodes[[src, dest], 1], zorder=3, c=cmap(i/n_routes))
            
        for i in range(self.n):
            plt.annotate(str(i), (self.nodes[i, 0], self.nodes[i, 1]))
        fig.patch.set_facecolor('white')
        plt.title(self.description['NAME'] + ' ' + self.current_formulation + ' Solution\nObjective = {:.2f}'.format(self.current_obj))
        plt.savefig('output/' + self.description['NAME'] + '_' + self.current_formulation, facecolor='white', dpi=100)
        
    def mtz(self):
        self.current_formulation = 'MTZ'
        model = Model('CVRP_MTZ')

        # Ranges
        I = [i for i in range(1, self.n)]
        I0 = [i for i in range(0, self.n)]

        # Variables
        x = {(i, j): model.binary_var('x_{0}_{1}'.format(i, j)) for i in I0
            for j in I0 if i!=j}
        u = {i: model.integer_var(0, self.Q, 'U_{0}'.format(i)) for i in I}

        # Objective
        model.minimize(model.sum(self.d[i, j] * x[(i, j)] for i in I0 for j in I0 if i!=j))

        # Subject to
        # (1)
        for i in I:
            model.add_constraint(model.sum(x[(j, i)] for j in I0 if i!=j) == 1)

        # (2)
        for i in I:
            model.add_constraint(model.sum(x[(i, j)] for j in I0 if i!=j) == 1)
        # (6)
        
        for i in I:
            for j in I:
                if i!=j:
                    model.add_constraint(u[i] - u[j] + self.Q * x[(i, j)] + (self.Q - self.q[i] - self.q[j]) * x[(j, i)] <= self.Q - self.q[j])

        # (7)

        for i in I:
            model.add_constraint(u[i] >= self.q[i] + model.sum(self.q[j] * x[(j, i)] for j in I if i!=j))
        
        # (8)

        for i in I:
            model.add_constraint(u[i] <= self.Q - model.sum(self.q[j] * x[(i, j)] for j in I if i!=j))

        # (9)

        for i in I:
            model.add_constraint(u[i] <= self.Q - (self.Q - self.q[i]) * x[(0, i)])
        
        # (10)
        
        for i in I:
            model.add_constraint(u[i] <= self.Q - (self.Q - np.delete(self.q, i).max()- self.q[i]) * x[(0, i)] - model.sum(self.q[j] * x[(i, j)] for j in I if j!=i))
        
        return model

    def gg(self):
        self.current_formulation = 'GG'
        model = Model('CVRP_GG')

        # Ranges
        I = [i for i in range(1, self.n)]
        I0 = [i for i in range(0, self.n)]
        # Variables
        x = {(i, j): model.binary_var('x_{0}_{1}'.format(i, j)) for i in I0
            for j in I0 if i!=j}

        F = {(i, j): model.continuous_var(0, name='F_{0}_{1}'.format(i, j)) for i in I0 for j in I0 if i!=j}

        # Objective

        model.minimize(model.sum(self.d[i, j] * x[(i, j)] for i in I0 for j in I0 if i!=j))

        # Subject to
        # (1)
        for i in I:
            model.add_constraint(model.sum(x[(j, i)] for j in I0 if i!=j) == 1)

        # (2)
        for i in I:
            model.add_constraint(model.sum(x[(i, j)] for j in I0 if i!=j) == 1)

        # (12)

        for i in I:
            model.add_constraint(model.sum(F[(j, i)] for j in I0 if i!=j) + self.q[i] == model.sum(F[(i, j)] for j in I0 if i!=j))
        
        # (13)

        for i in I0:
            for j in I0:
                if i!=j:
                    model.add_constraint(F[(i, j)] >= self.q[i] * x[(i, j)])
        
        # (14)

        for i in I0:
            for j in I0:
                if i!=j:
                    model.add_constraint(F[(i, j)] <= (self.Q - self.q[j]) * x[(i, j)])
        
        return model
    
    def bhm(self):
        _thisnodes = np.concatenate((self.nodes, self.nodes[:1, :]), axis=0)
        _d = np.rint(distance_matrix(_thisnodes, _thisnodes))
        self.current_formulation = 'BHM'
        model = Model('CVRP_BHM')

        # Ranges
        I = [i for i in range(1, self.n)]
        I0 = [i for i in range(0, self.n)]
        I0_p = [i for i in range(0, self.n + 1)]

        # Variables
        x = {(i, j): model.binary_var('x_{0}_{1}'.format(i, j)) for i in I0_p
            for j in I0_p}
        G = {(i, j): model.continuous_var(0, name='G_{0}_{1}'.format(i, j)) 
            for i in I0_p for j in I0_p}
        
        # FO
        model.minimize((model.sum(_d[i, j] * x[(i, j)] for i in I0 for j in I0_p if i!=j) + model.sum(_d[i, self.n] * x[(i, self.n)] for i in I))/2)
        # Subject to

        # (17)
        for i in I:
                model.add_constraint(model.sum((G[(j, i)] - G[(i, j)]) for j in I0_p) == 2 * self.q[i], ctname='R:17 (i={})'.format(i)) 

        # (18)
        model.add_constraint(model.sum(G[(0, j)] for j in I) == model.sum(self.q[j] for j in I), ctname='R: 18')

        # (19)

        model.add_constraint(model.sum(G[(j, 0)] for j in I) == self.K * self.Q - model.sum(self.q[j] for j in I), ctname='R:19')

        # (20)

        model.add_constraint(model.sum(G[(self.n, j)] for j in I) == self.K * self.Q, ctname='R20')

        # (21)

        for i in I0_p:
            for j in I0_p:
                if i != j and i != self.n:
                    model.add_constraint(G[(i, j)] + G[(j, i)] == self.Q * x[(i, j)], ctname='R:21 (i={0}, j={1})'.format(i, j))

        # (22)

        for i in I:
            model.add_constraint(model.sum(x[(i, j)] for j in I0_p if j > i) + model.sum(x[(i, j)] for j in I0_p if j < i) == 2, ctname='R:22 (i={0})'.format(i))
        
        return model
    
    def mcf(self):
        self.current_formulation = 'MCF'
        model = Model('CVRP_MCF')

        # Ranges
        I = [i for i in range(1, self.n)]
        I0 = [i for i in range(0, self.n)]
        # Variables
        x = {(i, j): model.binary_var('x_{0}_{1}'.format(i, j)) for i in I0
            for j in I0 if i!=j}
        
        F = {(i, j, k): model.binary_var('F_{0}_{1}_{2}'.format(i, j, k)) for i in I0 for j in I0 for k in I if i!=j}
        G = {(i, j, k): model.binary_var('G_{0}_{1}_{2}'.format(i, j, k)) for i in I0 for j in I0 for k in I if i!=j}

        # Objective

        model.minimize(model.sum(self.d[i, j] * x[(i, j)] for i in I0 for j in I0 if i!=j))

        # Subject to
        # (1)
        for i in I:
            model.add_constraint(model.sum(x[(j, i)] for j in I0 if i!=j) == 1)

        # (2)
        for i in I:
            model.add_constraint(model.sum(x[(i, j)] for j in I0 if i!=j) == 1)


        # (25)
        for k in I:
            model.add_constraint(model.sum(F[(0, j, k)] for j in I) - model.sum(F[(j, 0, k)] for j in I) == 1)
        
        # (26)

        for i in I:
            for k in I:
                if i!=k:
                    model.add_constraint(model.sum(F[i, j, k] for j in I0 if j!=i) - model.sum(F[j, i, k] for j in I0 if j!=i) == 0)

        # (27)

        for k in I:
            model.add_constraint(model.sum(G[(j, 0, k)] for j in I) - model.sum(G[(0, j, k)] for j in I) == 1)
        
        # (28)
        for i in I:
            for k in I:
                if i!=k:
                    model.add_constraint(model.sum(G[i, j, k] for j in I0 if j!=i) - model.sum(G[j, i, k] for j in I0 if j!=i) == 0)
        
        # (29)

        for i in I0:
            for j in I0:
                for k in I:
                    if i!=j:
                        model.add_constraint(F[(i, j, k)] + G[(i, j, k)] <= x[(i, j)])
        
        # (30)

        for i in I0:
            for j in I0:
                if i!=j:
                    model.add_constraint(model.sum(self.q[k] * (F[(i, j, k)] + G[(i, j, k)]) for k in I if (k!=i) and (k!=j)) <= (self.Q - self.q[i] - self.q[j]) * x[(i, j)])
        
        return model
        
    def solve(self, formulation='mtz', timelimit=120):
        formulations = {'mtz': self.mtz,
                        'gg': self.gg,
                        'bhm': self.bhm,
                        'mcf': self.mcf}

        model = formulations[formulation]()
        model.parameters.timelimit = timelimit
        sol = model.solve(log_output=True)



        self.current_obj = sol.objective_value
        
        if formulation == 'bhm':
            self.current_sol = {(i, j): model.get_var_by_name('x_{0}_{1}'.format(i, j)) for i in range(self.n + 1) for j in range(self.n + 1) if i!=j}
        else:
            self.current_sol = {(i, j): model.get_var_by_name('x_{0}_{1}'.format(i, j)) for i in range(self.n) for j in range(self.n) if i!=j}
        
        self.plot()

        sol_details = model.get_solve_details()
        output = open('output/' + self.description['NAME'] + '_' + formulation + '.txt', 'w')
        output.write('Optimal: ' + str(self.optimal) + '\n')
        output.write('Best Bound: ' + str(sol_details.best_bound) + '\n')
        output.write('Columns: ' + str(sol_details.columns) + '\n')
        output.write('N Iter: ' + str(sol_details.nb_iterations) + '\n')
        output.write('N Nodes: ' + str(sol_details.nb_nodes_processed) + '\n')
        output.write('Status: ' + str(sol_details.status) + '\n')
        output.write('Time: ' + str(sol_details.time) + '\n')
        output.write('Limit: ' + str(sol_details.has_hit_limit()) + '\n')
        output.write(str(sol))
        output.close()
        return sol


if __name__ == '__main__':

    set_p = ['http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n16-k8.vrp']
            #'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n19-k2.vrp',
            #'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n20-k2.vrp',
            #'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n21-k2.vrp',
            #'http://vrp.atd-lab.inf.puc-rio.br/media/com_vrp/instances/P/P-n22-k2.vrp']
    for url in set_p:
        instance = CVRP.from_url(url)
        for formulation in ['gg']:
            sol = instance.solve(formulation=formulation, timelimit=2000)