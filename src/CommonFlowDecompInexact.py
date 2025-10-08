import networkx as nx
import kCommonFlowDecompInexact as kCFDI

class CommonFlowDecompInexact:
    def __init__(self, G: nx.DiGraph, num_flows: int, maximum_k: int, flow_attr: str = "flow", subpath_constr: list = []):
        self.G = G
        self.num_flows = num_flows
        self.maximum_k = maximum_k
        self.flow_attr = flow_attr
        self.subpath_constr = subpath_constr
        
    def solve(self):
        for k in range(1,self.maximum_k):
            myDecomp = kCFDI.KCommonFlowDecompInexact(self.G, self.num_flows, k, self.flow_attr, self.subpath_constr)
            myDecomp.build_model()
            if myDecomp.solve_model():
                solution = f"Found a solution with {k} distinct paths:\n" + myDecomp.get_model_solution()
                return solution
        return "No solution found in specified range of k."