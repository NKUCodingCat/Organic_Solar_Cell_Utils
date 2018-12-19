
import mol2, networkx, re, argparse, collections

parser = argparse.ArgumentParser()
parser.add_argument('-i', help='Input, Mol2 file name', required=True)
parser.add_argument('-o', help='Output, GRO file name', required=True)

args = parser.parse_args()

def _(s):
    return re.match("^([^\.]+)(\..+)?$", s).group(1)

G=networkx.Graph()
D = mol2.read_Mol2_file(args.i)[0]

M = { i.num:(_(i.type), i.X/10.0, i.Y/10.0, i.Z/10.0, i.Q)  for i in D.atom_list}
N = [ (j.a1_num, j.a2_num) for j in D.bond_list]
G.add_nodes_from(M)
G.add_edges_from(N)

sub_graphs = list(networkx.connected_component_subgraphs(G))

print "find", len(sub_graphs), "molecules"
with open(args.o, "w") as f:
    print >> f , "I DONT KNOW WHAT I AM DOING..."
    print >> f , len(M)
    for i, sub_G in enumerate(sub_graphs):
        # Name = "%dUNL"%i
        print "RES-%5d "%(i+1), " ".join(["%s%d"%(k, v) for k, v in sorted(dict(collections.Counter([M[n][0] for n in sub_G.nodes])).items(), key=lambda x:x[0])])
        for n in sub_G.nodes:
            t = M[n]
            print  >> f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"%(i+1, "UNL", t[0], n, t[1], t[2], t[3])
    print >> f ,"   0.00000   0.00000   0.00000"