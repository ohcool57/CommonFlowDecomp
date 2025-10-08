import networkx as nx
import kCommonFlowDecompMinPathErr as kCFDPE

class CommonFlowDecompMinPathErr:
    def __init__(self, G: nx.DiGraph, num_flows: int, maximum_k: int, flow_attr: str = "flow", subpath_constr: list = []):
        self.G = G
        self.num_flows = num_flows
        self.maximum_k = maximum_k
        self.flow_attr = flow_attr
        self.subpath_constr = subpath_constr
        
    def solve(self):
        last_obj = float("inf")
        solution = ""
        for k in range(1,self.maximum_k+2):
            myDecomp = kCFDPE.KCommonFlowDecompMinPathErr(self.G, self.num_flows, k, self.flow_attr, self.subpath_constr)
            myDecomp.build_model()
            new_obj = myDecomp.solve_model()
            if new_obj == float("inf"):
                last_obj = new_obj
                solution = solution + f"No solution for {k} paths\n"
            elif new_obj < last_obj:
                last_obj = new_obj
                last_solution = myDecomp.get_model_solution()
                solution = solution + f"Found a solution with {k} distinct paths and total path error {last_obj}\n"
            elif new_obj == last_obj:
                solution = solution + f"Optimal solution: {k - 1} distinct paths and total path error {last_obj}:\n{last_solution}"
                return solution
        return "No solution found in specified range of k."