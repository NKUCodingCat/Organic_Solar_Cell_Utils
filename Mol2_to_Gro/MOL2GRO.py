
import mol2, networkx, re, argparse, collections, re

NAME_to_ELECTRON = {'BE': 4, 'BA': 56, 'BH': 107, 'BI': 83, 'BK': 97, 'BR': 35, 'RU': 44, 'RE': 75, 'RF': 104, 
    'RG': 111, 'RA': 88, 'RB': 37, 'RN': 86, 'RH': 45, 'TM': 69, 'H': 1, 'P': 15, 'GE': 32, 'GD': 64, 'GA': 31, 
    'UUT': 113, 'OS': 76, 'HS': 108, 'ZN': 30, 'HO': 67, 'HF': 72, 'HG': 80, 'HE': 2, 'PR': 59, 'PT': 78, 'PU': 94, 
    'UUO': 118, 'PB': 82, 'PA': 91, 'PD': 46, 'PO': 84, 'PM': 61, 'C': 6, 'K': 19, 'O': 8, 'S': 16, 'W': 74, 'EU': 63, 
    'ZR': 40, 'ER': 68, 'MD': 101, 'MG': 12, 'MO': 42, 'MN': 25, 'MT': 109, 'AU': 79, 'FE': 87, 'FL': 114, 'FM': 100, 
    'NI': 28, 'NO': 102, 'NA': 11, 'NB': 41, 'ND': 60, 'NE': 10, 'ES': 99, 'NP': 93, 'B': 5, 'CO': 27, 'CN': 112, 'CM': 96, 
    'CL': 17, 'CA': 20, 'CF': 98, 'CE': 58, 'CD': 48, 'V': 23, 'CS': 55, 'CR': 24, 'CU': 29, 'SR': 38, 'UUP': 115, 'UUS': 117, 
    'KR': 36, 'SI': 14, 'SN': 50, 'SM': 62, 'SC': 21, 'SB': 51, 'SG': 106, 'SE': 34, 'YB': 70, 'DB': 105, 'DY': 66, 'DS': 110, 
    'LA': 57, 'F': 9, 'LI': 3, 'LV': 116, 'TL': 81, 'LU': 71, 'LR': 103, 'TH': 90, 'TI': 22, 'TE': 52, 'TB': 65, 'TC': 43, 
    'TA': 73, 'AC': 89, 'AG': 47, 'I': 53, 'IR': 77, 'AM': 95, 'AL': 13, 'AS': 33, 'AR': 18, 'U': 92, 'AT': 85, 'IN': 49, 
    'Y': 39, 'N': 7, 'XE': 54}

GREP_ELETRON = lambda s: sum( map( lambda s1: (NAME_to_ELECTRON[s1[0]]*int(s1[1])), map(lambda s0: re.match("^([^\s\d]+)(\d+)$", s0).groups(), re.findall("[^\s\d]+\d+", s))) )

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

print "find", len(sub_graphs), "molecules in", 
# with open(args.o, "w") as f:
#     print >> f , "I DONT KNOW WHAT I AM DOING..."
#     print >> f , len(M)
#     for i, sub_G in enumerate(sub_graphs):
#         # Name = "%dUNL"%i
#         print "RES-%5d "%(i+1), " ".join(["%s%d"%(k, v) for k, v in sorted(dict(collections.Counter([M[n][0] for n in sub_G.nodes])).items(), key=lambda x:x[0])])
#         for n in sub_G.nodes:
#             t = M[n]
#             print  >> f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"%(i+1, "UNL", t[0], n, t[1], t[2], t[3])
#     print >> f ,"   0.00000   0.00000   0.00000"



DATA = [   # id, obj,     chemical formula
    (i+1, sub_G, " ".join(["%s%d"%(k, v) for k, v in sorted(dict(collections.Counter([M[n][0] for n in sub_G.nodes])).items(), key=lambda x:x[0])]))
    for i, sub_G in enumerate(sub_graphs)
]

Kind = set([j[2] for j in DATA])
print len(Kind), "kinds"
assert len(Kind) == 2, "Find MORE THAN 2 KINDS of molecule, please check"

Grouped_data = {i:list() for i in Kind}
for i in DATA: Grouped_data[i[2]].append(i)

with open(args.o, "w") as f:
    
    print >> f , "I DONT KNOW WHAT I AM DOING..."
    print >> f , len(M)

    CurAtom, CurResidue = 1, 1
    for i, j in enumerate(Grouped_data):
        NAME, mols = i, Grouped_data[j]
        print "Kind%2d (%5d e-) = %s"%(i, GREP_ELETRON(j), j)
        for mol in mols:
            for n in mol[1].nodes:
                t = M[n]
                print  >> f, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f"%(CurResidue, "MG%d"%NAME, t[0], CurAtom%10000, t[1], t[2], t[3])
                CurAtom += 1
            CurResidue  += 1
    print >> f ,"   0.00000   0.00000   0.00000"

############################
# GRO File FORMAT Ref: http://manual.gromacs.org/archive/5.0.3/online/gro.html

# Lines contain the following information (top to bottom):

#     - title string (free format string, optional time in ps after 't=')
#     - number of atoms (free format integer)
#     - one line for each atom (fixed format, see below)
#     - box vectors (free format, space separated reals), values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), 
#          the last 6 values may be omitted (they will be set to zero). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.

# This format is fixed, ie. all columns are in a fixed position. Optionally (for now only yet with trjconv) you can write gro 
# files with any number of decimal places, the format will then be n+5 positions with n decimal places (n+1 for velocities) instead of
# 8 with 3 (with 4 for velocities). Upon reading, the precision will be inferred from the distance between the decimal points (which will be n+5). 
# Columns contain the following information (from left to right):

#     - residue number (5 positions, integer)
#     - residue name (5 characters)
#     - atom name (5 characters)
#     - atom number (5 positions, integer)
#     - position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
#     - velocity (in nm/ps (or km/s), x y z in 3 columns, each 8 positions with 4 decimal places)

# Note that separate molecules or ions (e.g. water or Cl-) are regarded as residues. If you want to write such a file in your 
# own program without using the GROMACS libraries you can use the following formats:

# C format   "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"

# Fortran format  (i5,2a5,i5,3f8.3,3f8.4)

# Note that this is the format for writing, as in the above example fields may be written without spaces, and therefore can not be read with the same format statement in C.