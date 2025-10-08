import networkx as nx
import gurobipy as gb
import utils

class KCommonFlowDecompBoundedErr:
    def __init__(self, G: nx.DiGraph, num_flows: int, k: int, error_bound: float, flow_attr: str = "flow", subpath_constr: list = []):
        if not nx.is_directed_acyclic_graph(G):
            print("uh oh")
            raise ValueError('Input graph is not a directed acyclic graph')
        if not utils.check_st_graph(G):
            print("uh oh")
            raise ValueError('Input graph is not an st graph')
        if not utils.check_correct_num_flows(G, num_flows, flow_attr):
            print("uh oh")
            raise ValueError('Number of flows does not match')
        if not utils.check_valid_flow_format(G, num_flows, flow_attr):
            print("uh oh")
            raise ValueError('Flow value must be int or float')
        if subpath_constr and not utils.check_subpath_constr(G, subpath_constr):
            print("uh oh")
        self.model = gb.Model()
        self.G = G
        self.num_flows = num_flows
        self.k = k
        self.flow_attr = flow_attr
        self.error_bound = error_bound
        self.w_max = utils.get_max_flow(self.G, self.num_flows, self.flow_attr)
        self.subpath_constr = [[] for _ in range(len(subpath_constr))]
        for subpath in range(len(subpath_constr)):
            for node in range(1,len(subpath_constr[subpath])):
                self.subpath_constr[subpath].append((subpath_constr[subpath][node-1], subpath_constr[subpath][node]))

        self.path_indexes = [(i, j) for i in range(self.k) for j in range(self.num_flows)]
        self.edge_indexes = [(u, v, i) for u, v in self.G.edges() for i in range(self.k)]
        self.edge_flows = {(u, v, j): data[self.flow_attr][j] for u, v, data in self.G.edges(data=True) for j in
                           range(self.num_flows)}
        self.pi_indexes = [(u, v, i, j) for u, v in self.G.edges() for i in range(self.k) for j in
                           range(self.num_flows)]
        self.subpath_indexes = [(i, p) for i in range(self.k) for p in range(len(self.subpath_constr))]

    def build_model(self):

        self.variable_name_prefixes = []

        self.path_vars = self.add_variables(indexes=self.path_indexes, name_prefix='w', ub=self.w_max)
        self.edge_vars = self.add_variables(indexes=self.edge_indexes, name_prefix='x', var_type="binary")
        self.pi_vars = self.add_variables(indexes=self.pi_indexes, name_prefix='pi', ub=self.w_max)
        self.subpath_vars = self.add_variables(indexes=self.subpath_indexes, name_prefix='r', var_type="binary")

        for v in self.G.nodes():
            predecessors = list(self.G.predecessors(v))
            successors = list(self.G.neighbors(v))
            if len(predecessors) == 0:
                for i in range(self.k):
                    self.model.addConstr(gb.quicksum(self.edge_vars[v, w, i] for w in successors) == 1,
                                         name=f"single_path_i={i}")
            elif len(successors) != 0:
                for i in range(self.k):
                    self.model.addConstr(gb.quicksum(self.edge_vars[u, v, i] for u in predecessors) ==
                                         gb.quicksum(self.edge_vars[v, w, i] for w in successors),
                                         name=f"flow_cons_v={v}_i={i}")
        for u, v in self.G.edges():
            for j in range(self.num_flows):
                self.model.addConstr(gb.quicksum(self.pi_vars[u, v, i, j] for i in range(self.k)) <=
                                     self.edge_flows[u, v, j] + self.error_bound,
                                     name=f"correct_flow_u={u}_v={v}_j={j}")
                self.model.addConstr(gb.quicksum(self.pi_vars[u, v, i, j] for i in range(self.k)) >=
                                     self.edge_flows[u, v, j] - self.error_bound,
                                     name=f"correct_flow_u={u}_v={v}_j={j}")

                for i in range(self.k):
                    self.add_binary_continuous_product_constraint(binary_var=self.edge_vars[u, v, i],
                                                                  continuous_var=self.path_vars[i, j],
                                                                  product_var=self.pi_vars[u, v, i, j], lb=0,
                                                                  ub=self.w_max, name=f"pi_u={u}_v={v}_i={i}_j={j}")



        ###PRIMARY FORMULATION -- EACH SUBPATH CONSTRAINT SATISFIED BY A SINGLE FLOW
        if self.subpath_constr:
            for p in range(len(self.subpath_constr)):
                self.model.addConstr(gb.quicksum(self.subpath_vars[i,p] for i in range(self.k)) >= 1,
                                     name=f"subpath_claim_p={p}")
                for i in range(self.k):
                    self.model.addConstr(gb.quicksum(self.edge_vars[u,v,i] for u, v in self.subpath_constr[p]) >=
                                         len(self.subpath_constr[p]) * self.subpath_vars[i,p],
                                         name=f"subpath_proof_i={i}_p={p}")

            for i in range(self.k):
                self.model.addConstr(gb.quicksum(self.path_vars[i,j] for j in range(self.num_flows)) >= 1,
                                     name=f"path_used_i={i}")

        ###ALTERNATIVE FORMULATION -- EACH SUBPATH CONSTRAINT SATISFIED BY ALL FLOWS
        # if self.subpath_constr:
        #     for j in range(self.num_flows):
        #         for p in range(len(self.subpath_constr)):
        #             self.model.addConstr(gb.quicksum(self.subpath_vars[i,j,p] for i in range(self.k)) >= 1,
        #                                  name=f"subpath_flow_claim_j={j}_p={p}")
        #             for i in range(self.k):
        #                 self.model.addConstr(gb.quicksum(self.edge_vars[u,v,i] for u, v in self.subpath_constr[p]) >=
        #                                      (len(self.subpath_constr[p]) - 1) * self.subpath_vars[i,j,p],
        #                                      name=f"subpath_proof_i={i}_j={j}_p={p}")
        #
        #         for i in range(self.k):
        #             self.model.addConstr(self.path_vars[i,j] >= self.path_vars[i,j,p],
        #                                  name=f"path_flow_used_i={i}_j={j}_p={p}")

    def solve_model(self):
        self.model.optimize()
        if self.model.status == gb.GRB.Status.OPTIMAL:
            return True
        else:
            return False

    def get_model_solution(self):
        solution = ""
        for i in range(self.k):
            solution = solution + f"Path {i+1} (carries weight"
            for j in range(self.num_flows):
                if j == self.num_flows - 1:
                    solution = solution + " and"
                solution = solution + f" {self.path_vars[(i, j)].X} for flow {j+1}"
                if j < self.num_flows - 1 and self.num_flows > 2:
                    solution = solution + ","
            solution = solution + "):\n"
            for u,v in self.G.edges():
                if self.edge_vars[u, v, i].X != 0:
                    solution = solution + f"{u}, "
            solution = solution + "t\n"
        for i in range(self.k):
            for p in range(len(self.subpath_constr)):
                if self.subpath_vars[(i,p)].X != 0:
                    solution = solution + f"Path {i+1} satisfies constraint {j}\n"
        return solution


    def add_variables(self, indexes, name_prefix: str, lb=0, ub=1, var_type="integer"):
        for prefix in self.variable_name_prefixes:
            if prefix.startswith(name_prefix) or name_prefix.startswith(prefix):
                print("uh oh")
                raise ValueError(
                    f"Variable name prefix {name_prefix} conflicts with existing variable name prefix {prefix}. "
                    f"Use a different name prefix."
                )

        self.variable_name_prefixes.append(name_prefix)

        var_type_map = {
            "integer": gb.GRB.INTEGER,
            "continuous": gb.GRB.CONTINUOUS,
            "binary": gb.GRB.BINARY,
        }
        vars = {}
        for index in indexes:
            vars[index] = self.model.addVar(
                lb=lb,
                ub=ub,
                vtype=var_type_map[var_type],
                name=f"{name_prefix}{index}",
            )
        self.model.update()
        return vars

    def add_binary_continuous_product_constraint(self, binary_var, continuous_var, product_var, lb, ub, name: str):
        self.model.addConstr(product_var <= ub * binary_var, name=name + "_a")
        self.model.addConstr(product_var >= lb * binary_var, name=name + "_b")
        self.model.addConstr(product_var <= continuous_var - lb * (1 - binary_var), name=name + "_c")
        self.model.addConstr(product_var >= continuous_var - ub * (1 - binary_var), name=name + "_d")