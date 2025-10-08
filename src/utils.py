import networkx as nx

def check_valid_multi_flow(G: nx.Graph, num_flows: int, perfect_flow: str = "perfect", flow_attr: str = "flow", subpath_constr: list = []) -> bool:
    if not nx.is_directed_acyclic_graph(G):
        print("uh oh")
        raise ValueError('Input graph is not a directed acyclic graph')
        return False
    if not check_st_graph(G):
        print("uh oh")
        raise ValueError('Input graph is not an st graph')
        return False
    if not check_correct_num_flows(G, num_flows, flow_attr):
        print("uh oh")
        raise ValueError('Number of flows does not match')
        return False
    if perfect_flow and not check_multi_flow_conservation(G, num_flows, flow_attr):
        print("uh oh")
        raise ValueError('Input graph does not conserve flow')
        return False
    if subpath_constr and not check_subpath_constr(G, subpath_constr):
        print("uh oh")
        return False
    return True

def check_st_graph(G: nx.DiGraph) -> bool:
    if not nx.is_directed_acyclic_graph(G): return False
    single_source = False
    single_sink = False
    for v in G.nodes():
        if G.in_degree(v) == 0:
            if single_source:
                return False
            else:
                single_source = True
                continue
        if G.out_degree(v) == 0:
            if single_sink:
                return False
            else:
                single_sink = True
                continue
    return single_source and single_sink

def check_correct_num_flows(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> bool:
    for u, v, data in G.edges(data = True):
        if len(data.get(flow_attr)) != num_flows:
            return False
    return True

def check_valid_flow_format(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> bool:
    for u, v, data in G.edges(data = True):
        if not all(isinstance(flow_val, (int, float)) for flow_val in data.get(flow_attr)):
            return False
    return True

def check_multi_flow_conservation(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> bool:
    for v in G.nodes():
        if G.out_degree(v) == 0 or G.in_degree(v) == 0:
            continue
        for j in range(num_flows):
            out_flow = 0
            for x, y, data in G.out_edges(v, data=True):
                out_flow += data[flow_attr][j]

            in_flow = 0
            for x, y, data in G.in_edges(v, data=True):
                in_flow += data[flow_attr][j]

            if out_flow != in_flow:
                return False
    return True

def check_subpath_constr(G: nx.DiGraph, subpath_constr: list) -> bool:
    if not isinstance(subpath_constr, list) or not all(
            isinstance(subpath, list) and all(isinstance(item, str) for item in subpath)
            for subpath in subpath_constr
    ):
        raise ValueError("data must be a list of lists of strings")
        return False

    for subpath in subpath_constr:
        for i in range(1,len(subpath)):
            if not G.has_edge(subpath[i-1], subpath[i]):
                raise ValueError("subpaths must be connected")
                return False

    return True

def check_valid_inexact_flows(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> bool:
    for u,v,data in G.edges(data=True):
        for j in range(num_flows):
            if not isinstance(data.get(flow_attr)[j], tuple):
                raise ValueError("Flow attributes must be bounds expressed as tuples")
                return False
            if len(data.get(flow_attr)[j]) != 2:
                raise ValueError("There must be a single upper and lower bound for each edge flow value")
                return False
            if not data.get(flow_attr)[j][0] <= data.get(flow_attr)[j][1]:
                raise ValueError("Lower bound must be less than or equal to upper bound for each edge flow value")
                return False
    return True

def get_max_flow(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> int:
    w_max = float("-inf")
    if not check_correct_num_flows(G, num_flows, flow_attr):
        print("uh oh")
        raise ValueError(
            "Some edges missing flows"
        )
    for u, v, data in G.edges(data=True):
        if not flow_attr in data:
            print("uh oh")
            raise ValueError(
                f"Edge ({u},{v}) does not have the required flow attribute '{flow_attr}'. Check that the attribute passed under 'flow_attr' is present in the edge data."
            )
        if any(flow < 0 for flow in data[flow_attr]):
            print("uh oh")
            raise ValueError(
                f"Edge ({u},{v}) has negative flow value {data[flow_attr]}. All flow values must be >=0."
            )
        w_max = max(w_max, max(data[flow_attr]))
    return w_max

def get_max_inexact_flow(G: nx.DiGraph, num_flows: int, flow_attr: str = "flow") -> int:
    w_max = float("-inf")
    if not check_correct_num_flows(G, num_flows, flow_attr):
        print("uh oh")
        raise ValueError(
            "Some edges missing flows"
        )
    for u, v, data in G.edges(data=True):
        if not flow_attr in data:
            print("uh oh")
            raise ValueError(
                f"Edge ({u},{v}) does not have the required flow attribute '{flow_attr}'. Check that the attribute passed under 'flow_attr' is present in the edge data."
            )
        if any(flow[0] < 0 for flow in data[flow_attr]):
            print("uh oh")
            raise ValueError(
                f"Edge ({u},{v}) has negative flow value {data[flow_attr]}. All flow values must be >=0."
            )
        w_max = max(w_max, max(data[flow_attr][1]))
    return w_max