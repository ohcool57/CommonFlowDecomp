import networkx as nx
import gurobipy as gb
import utils

class KCommonFlowDecompMinPathErr:
    def __init__(self, G: nx.DiGraph, num_flows: int, k: int, flow_attr: str = "flow", subpath_constr: list = []):
        if not nx.is_directed_acyclic_graph(G):
            print("uh oh")
            raise ValueError('Input graph is not a directed acyclic graph')
        if not utils.check_st_graph(G):
            print("uh oh")
            raise ValueError('Input graph is not an st graph')
        if not utils.check_valid_flow_format(G, num_flows, flow_attr):
            print("uh oh")
            raise ValueError('Flow value must be int or float')
        if not utils.check_correct_num_flows(G, num_flows, flow_attr):
            print("uh oh")
            raise ValueError('Number of flows does not match')
        if subpath_constr and not utils.check_subpath_constr(G, subpath_constr):
            print("uh oh")
        self.model = gb.Model()
        self.G = G
        self.num_flows = num_flows
        self.k = k
        self.flow_attr = flow_attr
        self.w_max = utils.get_max_flow(self.G, self.num_flows, self.flow_attr)

        self.path_indexes = [(i, j) for i in range(self.k) for j in range(self.num_flows)]
        self.edge_indexes = [(u, v, i) for u, v in self.G.edges() for i in range(self.k)]
        self.edge_flows = {(u, v, j): data[self.flow_attr][j] for u, v, data in self.G.edges(data=True) for j in
                           range(self.num_flows)}
        self.pi_indexes = [(u, v, i, j) for u, v in self.G.edges() for i in range(self.k) for j in
                           range(self.num_flows)]


    def build_model(self):

        self.variable_name_prefixes = []

        self.path_vars = self.add_variables(indexes=self.path_indexes, name_prefix='w', ub=self.w_max)
        # self.path_slack_vars = self.add_variables(indexes=self.path_indexes, name_prefix='rho', ub=self.w_max)
        # self.gamma_vars = self.add_variables(indexes=self.pi_indexes, name_prefix='gamma', ub=self.w_max)
        # USE THE ABOVE TWO LINES IF USING A PATH-FLOW SLACK. USE THE BELOW TWO LINES IF USING A PATH SLACK.
        self.path_slack_vars = self.add_variables(indexes=range(self.k), name_prefix='rho', ub=self.w_max)
        self.gamma_vars = self.add_variables(indexes=self.edge_indexes, name_prefix='gamma', ub=self.w_max)
        self.edge_vars = self.add_variables(indexes=self.edge_indexes, name_prefix='x', var_type="binary")
        self.pi_vars = self.add_variables(indexes=self.pi_indexes, name_prefix='pi', ub=self.w_max)

        for v in self.G.nodes():
            predecessors = list(self.G.predecessors(v))
            successors = list(self.G.neighbors(v))
            if len(predecessors) == 0:
                for i in range(self.k):
                    self.model.addConstr(gb.quicksum(self.edge_vars[v, w, i] for w in successors) == 1, name=f"single_path_i={i}")
            elif len(successors) != 0:
                for i in range(self.k):
                    self.model.addConstr(gb.quicksum(self.edge_vars[u, v, i] for u in predecessors) == gb.quicksum(self.edge_vars[v, w, i] for w in successors), name=f"flow_cons_v={v}_i={i}")

        for u, v in self.G.edges():
            for j in range(self.num_flows):
                for i in range(self.k):
                    self.add_binary_continuous_product_constraint(binary_var=self.edge_vars[u, v, i], continuous_var=self.path_vars[i, j], product_var=self.pi_vars[u, v, i, j], lb=0, ub=self.w_max, name=f"pi_u={u}_v={v}_i={i}_j={j}")
                    # self.add_binary_continuous_product_constraint(binary_var=self.edge_vars[u, v, i], continuous_var=self.path_slack_vars[i, j], product_var=self.gamma_vars[u, v, i, j], lb=0, ub=self.w_max, name=f"gamma_u={u}_v={v}_i={i}_j={j}")
                    # USE THE ABOVE LINE IF USING A PATH-FLOW SLACK. USE THE BELOW LINE IF USING A PATH SLACK.
                    self.add_binary_continuous_product_constraint(binary_var=self.edge_vars[u, v, i], continuous_var=self.path_slack_vars[i], product_var=self.gamma_vars[u, v, i], lb=0, ub=self.w_max, name=f"gamma_u={u}_v={v}_i={i}_j={j}")
                # self.model.addConstr(self.edge_flows[u,v,j] - gb.quicksum(self.pi_vars[u,v,i,j] for i in range(self.k)) <= gb.quicksum(self.gamma_vars[u,v,i,j] for i in range(self.k)), name=f"path_slack_a_u={u}_v={v}_j={j}")
                # self.model.addConstr(self.edge_flows[u,v,j] - gb.quicksum(self.pi_vars[u,v,i,j] for i in range(self.k)) >= - gb.quicksum(self.gamma_vars[u,v,i,j] for i in range(self.k)), name=f"path_slack_b_u={u}_v={v}_j={j}")
                # USE THE ABOVE TWO LINES IF USING A PATH-FLOW SLACK. USE THE BELOW TWO LINES IF USING A PATH SLACK.
                self.model.addConstr(self.edge_flows[u,v,j] - gb.quicksum(self.pi_vars[u,v,i,j] for i in range(self.k)) <= gb.quicksum(self.gamma_vars[u,v,i] for i in range(self.k)), name=f"path_slack_a_u={u}_v={v}_j={j}")
                self.model.addConstr(self.edge_flows[u,v,j] - gb.quicksum(self.pi_vars[u,v,i,j] for i in range(self.k)) >= - gb.quicksum(self.gamma_vars[u,v,i] for i in range(self.k)), name=f"path_slack_b_u={u}_v={v}_j={j}")
        # self.model.setObjective(gb.quicksum(self.path_slack_vars[i,j] for i in range(self.k) for j in range(self.num_flows)))
        # USE THE ABOVE LINE IF USING A PATH-FLOW SLACK. USE THE BELOW LINE IF USING A PATH SLACK.
        self.model.setObjective(gb.quicksum(self.path_slack_vars[i] for i in range(self.k)))

    def solve_model(self):
        if self.model.status == gb.GRB.Status.OPTIMAL:
            return self.model.ObjVal
        else:
            self.model.optimize()
            if self.model.status == gb.GRB.Status.OPTIMAL:
                return self.model.ObjVal
            else:
                return float("inf")

    def get_model_solution(self):
        solution = ""
        for i in range(self.k):
            solution = solution + f"Path {i+1} (carries weight "
            for j in range(self.num_flows):
                if j == self.num_flows - 1:
                    solution = solution + "and "
                # solution = solution + f"carries weight {self.path_vars[(i, j)].X} and slack {self.path_slack_vars[(i,j)].X} for flow {j+1}"
                # USE THE ABOVE LINE IF USING A PATH-FLOW SLACK. USE THE BELOW LINE IF USING A PATH SLACK.
                solution = solution + f"{self.path_vars[(i, j)].X} for flow {j+1}"
                if j < self.num_flows - 1 and self.num_flows > 2:
                    solution = solution + ", "
            # solution = solution + "):\n"
            # USE THE ABOVE LINE IF USING A PATH-FLOW SLACK. USE THE BELOW LINE IF USING A PATH SLACK.
            solution = solution + f", with slack {self.path_slack_vars[i].X}):\n"
            for u,v in self.G.edges():
                if self.edge_vars[u, v, i].X != 0:
                    solution = solution + f"{u}, "
            solution = solution + "t\n"
        return solution


    def add_variables(self, indexes, name_prefix: str, lb=0, ub=1, var_type="integer"):
        for prefix in self.variable_name_prefixes:
            if prefix.startswith(name_prefix) or name_prefix.startswith(prefix):
                print("uh oh")
                raise ValueError(
                    f"Variable name prefix {name_prefix} conflicts with existing variable name prefix {prefix}. Use a different name prefix."
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