import networkx as nx
import kCommonFlowDecompMinErr as kCFDME

class CommonFlowDecompMinErr:
    def __init__(self, G: nx.DiGraph, num_flows: int, maximum_k: int, flow_attr: str = "flow", subpath_constr: list = []):
        self.G = G
        self.num_flows = num_flows
        self.maximum_k = maximum_k
        self.flow_attr = flow_attr
        self.subpath_constr = subpath_constr
        
    def solve(self, output: bool = False):
        last_obj = float("inf")
        for k in range(1,self.maximum_k+2):
            myDecomp = kCFDME.KCommonFlowDecompMinErr(self.G, self.num_flows, k, self.flow_attr, self.subpath_constr)
            myDecomp.build_model()
            new_obj = myDecomp.solve_model()
            paths = myDecomp.get_model_paths()
            del myDecomp
            if new_obj < last_obj:
                last_obj = new_obj
            elif new_obj == last_obj:
                if output:
                    print(f"Optimal solution: {k - 1} distinct paths and total error {last_obj}:\n{last_solution}")
                return paths
            if new_obj == 0:
                if output:
                    print(f"Optimal solution: {k} distinct paths and total error {last_obj}:\n{last_solution}")
                return paths
        if output:
            print("No optimal solution found in specified range of k.")
        return paths