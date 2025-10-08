import networkx as nx
import kCommonFlowDecompBoundedErr as kCFDBE


class CommonFlowDecompBoundedErr:
    def __init__(self, G: nx.DiGraph, num_flows: int, maximum_k: int, error_bound: float, flow_attr: str = "flow",
                 subpath_constr: list = []):
        self.G = G
        self.num_flows = num_flows
        self.maximum_k = maximum_k
        self.error_bound = error_bound
        self.flow_attr = flow_attr
        self.subpath_constr = subpath_constr

    def solve(self):
        for k in range(1, self.maximum_k):
            myDecomp = kCFDBE.KCommonFlowDecompBoundedErr(G=self.G, num_flows = self.num_flows, k=k, error_bound = self.error_bound, flow_attr=self.flow_attr, subpath_constr=self.subpath_constr)
            myDecomp.build_model()
            if myDecomp.solve_model():
                solution = f"Found a solution with {k} distinct paths:\n" + myDecomp.get_model_solution()
                return solution
        return "No solution found in specified range of k."