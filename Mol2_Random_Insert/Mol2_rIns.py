import numpy, argparse, networkx, copy, random
import mol2


# As = [(i.num, i.name, i.X, i.Y, i.Z) for i in Q.atom_list]
# Es = [(i.a1_num, i.a2_num, i.type)   for i in Q.bond_list]

# print As

# G = networkx.Graph()
# for i, v in enumerate(As):
#     G.add_node(v[0], p = v)

# for i, v in enumerate(Es):
#     G.add_edge(v[0], v[1], attr_dict={"typ": v[2]})

# print [ list(n[1]["p"][2:]) for n in G.nodes(data=True)]

class mol2_cls(object):
    
    def __init__(self, mol2obj):

        mol2obj = copy.deepcopy(mol2obj)
        As = [(i.num, i, i.X, i.Y, i.Z)    for i in mol2obj.atom_list]
        Es = [(i.a1_num, i.a2_num, i.type) for i in mol2obj.bond_list]
        G = networkx.Graph()

    
        
        self.Cord = numpy.array( [i[2:] for i in As] )
        Coord_avg = self.Cord.mean(axis = 0)
        self.Cord = self.Cord - Coord_avg # Centralize

        for i, v in enumerate(As):
            G.add_node(v[0], p = (v[0], v[1], v[2]-Coord_avg[0], v[3]-Coord_avg[1], v[4]-Coord_avg[2]))

        for i, v in enumerate(Es):
            G.add_edge(v[0], v[1], attr_dict={"typ": v[2]})
        
        self.Grph = G
        self.Atom_Count = len(As)
        self.Tran = numpy.array([0., 0., 0.])
        self.Rota = numpy.array([0., 0., 0.])
        self.Atom_Radii   = numpy.array([self._Atom_to_Radii(a[1].name) for a in As])
        self.Farthest_Dis = max( ((self.Cord**2).sum(axis = 1))**.5 + self.Atom_Radii )

        # self.Rad

    def Get_Cur_Cord(self):

        v  = self._Rot(self.Cord, self.Rota[0], self.Rota[1], self.Rota[2])
        v += self.Tran
        return v

    def Get_Cur_Center(self):
        return self.Tran 

    def Mol_Grph_w_Absolute_Cord(self):
        Ns = self.Grph.nodes(data=True)
        Es = self.Grph.edges(data=True)

        L = networkx.Graph()
        for i in Ns:
            L.add_node(i[0], C = (i[1]["p"][1], self._Rot(numpy.array(i[1]["p"][2:]), self.Rota[0], self.Rota[1], self.Rota[2]) + self.Tran))
        for j in Es:
            L.add_edge(j[0], j[1], attr_dict=j[2]["attr_dict"])
        return L

    @staticmethod
    def _Rot(x, Rx, Ry, Rz):
        Sx, Cx = numpy.sin(Rx), numpy.cos(Rx)
        Sy, Cy = numpy.sin(Ry), numpy.cos(Ry)
        Sz, Cz = numpy.sin(Rz), numpy.cos(Rz)
        _      = numpy.dot
        

        R_m = numpy.array(
            [
                [     Cy*Cz,          -Cy*Sz,          Sy ],
                [Cz*Sx*Sy+Cx*Sz,   Cx*Cz-Sx*Sy*Sz,  -Cy*Sx],
                [-Cx*Cz*Sy+Sx*Sz,  Cz*Sx+Cx*Sy*Sz,   Cx*Cy],
            ]
        )

        # Geo_Center = x.mean(axis=0)
        # Tmp = x - Geo_Center
        # Tmp = _(Tmp, R_m) + Geo_Center 
        return _(x, R_m)

    @staticmethod
    def _Atom_to_Radii(Atom_name):
        T = {"H" : 0.38, "He": 0.32, "Li": 1.34, "Be": 0.90, "B" : 0.82, "C" : 0.77, "N" : 0.75, "O" : 0.73, "F" : 0.71, "Ne": 0.69, "Na": 1.54, \
        "Mg": 1.30, "Al": 1.18, "Si": 1.11, "P" : 1.06, "S" : 1.02, "Cl": 0.99, "Ar": 0.97, "K" : 1.96, "Ca": 1.74, "Sc": 1.44, "Ti": 1.36, "V" : 1.25, \
        "Cr": 1.27, "Mn": 1.39, "Fe": 1.25, "Co": 1.26, "Ni": 1.21, "Cu": 1.38, "Zn": 1.31, "Ga": 1.26, "Ge": 1.22, "As": 1.19, "Se": 1.16, "Br": 1.14, \
        "Kr": 1.10, "Rb": 2.11, "Sr": 1.92, "Y" : 1.62, "Zr": 1.48, "Nb": 1.37, "Mo": 1.45, "Tc": 1.56, "Ru": 1.26, "Rh": 1.35, "Pd": 1.31, "Ag": 1.53, \
        "Cd": 1.48, "In": 1.44, "Sn": 1.41, "Sb": 1.38, "Te": 1.35, "I" : 1.33, "Xe": 1.30, "Cs": 2.25, "Ba": 1.98, "La": 1.69, "Ce": 1.818,"Pr": 1.824,\
        "Nd": 1.814,"Pm": 1.834,"Sm": 1.804,"Eu": 1.804,"Gd": 1.804,"Tb": 1.773,"Dy": 1.781,"Ho": 1.762,"Er": 1.761,"Tm": 1.759,"Yb": 1.76, "Lu": 1.60, \
        "Hf": 1.50, "Ta": 1.38, "W" : 1.46, "Re": 1.59, "Os": 1.28, "Ir": 1.37, "Pt": 1.28, "Au": 1.44, "Hg": 1.49, "Tl": 1.48, "Pb": 1.47, "Bi": 1.46}

        try:
            return T[Atom_name]
        except:
            print "Error Occured when we processing Atom name \"%s\""%Atom_name
            raise

class Box(object):

    delta = numpy.array([(i, j, k) for i in [-1, 0, 1] for j in [-1, 0, 1] for k in [-1, 0, 1] ])

    def __init__(self, size):
        self.mol_list = []
        size = numpy.array(size).flatten()
        assert size.size == 3
        self.size = size

    def Add_mol(self, mol_obj):
        if not self.Check_Collision(mol_obj):
            return False
        self.mol_list.append(copy.deepcopy(mol_obj))
        return True
    
    def Check_Collision(self, new_obj):
        # TODO: add boundary support
        __OLD = copy.deepcopy(new_obj.Tran)

        for d in self.delta:

            new_obj.Tran = copy.deepcopy(__OLD) + d*self.size
            for mol_F in self.mol_list:
                if not self.Coarse_Check(mol_F, new_obj):
                    continue
                if not self.Fine_Check(mol_F, new_obj):
                    continue
                return False

        new_obj.Tran = __OLD
        return True

    @staticmethod   # when collision occured , return True
    def Coarse_Check(mol_Fix, mol_Ins):
        if ((mol_Ins.Get_Cur_Center() - mol_Fix.Get_Cur_Center())**2).sum() < ( mol_Ins.Farthest_Dis + mol_Fix.Farthest_Dis )**2 :
           return True
        return False

    @staticmethod   # when collision occured , return True
    def Fine_Check(mol_Fix, mol_Ins):

        _C1 = mol_Fix.Get_Cur_Cord()
        _C2 = mol_Ins.Get_Cur_Cord()
        _A1 = mol_Fix.Atom_Radii
        _A2 = mol_Ins.Atom_Radii  # DAMN IT!

        _I1 = numpy.array(range(_C1.shape[0]))
        _I2 = numpy.array(range(_C2.shape[0]))
        
        _Ind1, _Ind2 = map(lambda x: x.flatten(), numpy.meshgrid(_I1, _I2))
        # print _C1[_I1].shape, _C2[_I2].shape, ( _A1[_I1] + _A2[_I2] ).shape
        # print _A1[_I1] + _A2[_I2]
        Q = (((_C1[_Ind1] - _C2[_Ind2])**2).sum(axis=1)) ** 0.5
        # ==========
        # __A__ = Q.argsort()[:20]
        # print Q[__A__], _Ind1[__A__], _Ind2[__A__]
        # print ( _A1[_Ind1] + _A2[_Ind2] )[Q.argsort()[:20]]
        # ==========
        Q -= ( _A1[_Ind1] + _A2[_Ind2] )
        return numpy.any( Q < .05 ) # for safety

    def Dump_box(self):

        ABSOLUTE_INDEX = 1
        BOX_ATOMS = []  #(atom_glo_index, atom_obj, atom_X, atom_Y, atom_Z)
        BOX_BONDS = []  #(a1_glo_index, a2_glo_index, type)

        assert len(set(map(id, self.mol_list))) == len(self.mol_list)

        for RES_i, mol in enumerate(self.mol_list):
            atom_ind_map = {}
            Q = mol.Mol_Grph_w_Absolute_Cord()
            for ind_atom, atom_in_local in enumerate(Q.nodes(data=True)):
                
                atom_ind_map[atom_in_local[0]] = ABSOLUTE_INDEX
                atom_prop, coord = atom_in_local[1]["C"]
                atom_prop.X,atom_prop.Y,atom_prop.Z = coord.tolist()
                atom_prop.num = ABSOLUTE_INDEX
                atom_prop.resnum = RES_i
                atom_prop.resname = "UN%d"%RES_i

                BOX_ATOMS.append(atom_prop)
                ABSOLUTE_INDEX += 1

            for b in Q.edges(data=True):

                BOX_BONDS.append(
                    (atom_ind_map[b[0]], atom_ind_map[b[1]], b[2]["attr_dict"]["typ"])
                )

        return BOX_ATOMS, BOX_BONDS, RES_i + 1, self.size

    def Box_2_mol_String(self, Mol_name): 
        A, B, Res, Box_size = self.Dump_box()
        Res_S = ""
        Res_S +=        "# A Comment line, Generated by HOME-MADE random insert. Developed by NKUCodingCat"
        Res_S += "\n" + "@<TRIPOS>MOLECULE"
        Res_S += "\n" + "%s"%Mol_name
        Res_S += "\n" + "%-5d %-5d %-5d 0     0"%(len(A), len(B), Res)
        Res_S += "\n" + "SMALL"
        Res_S += "\n" + "USER_CHARGES\n\n"
        Res_S += "\n" + "@<TRIPOS>ATOM"
        for i in A:
            Res_S += "\n" + "%-6d %-4s %9.4f %9.4f %9.4f %-5s %4s %6s %9.4f" % (
                i.num, i.name, i.X, i.Y, i.Z, i.type, i.resnum, i.resname, i.Q
            )
        Res_S += "\n" + "@<TRIPOS>BOND"
        for j, v in enumerate(B):
            Res_S += "\n" + "%-5d %-5d %-5d %s"%((j+1, )+v)
        Res_S += "\n" + "@<TRIPOS>CRYSIN"
        Res_S += "\n" + "%-12.8f  %-12.8f  %-12.8f  90   90   90   1   1"%tuple(Box_size.tolist())

        return Res_S

# TODO: we need argparser, maybe some optimization for checking collision

pc71 =  mol2_cls(mol2.read_Mol2_file("pc71bm.mol2")[0]         ) # 266
drcn =  mol2_cls(mol2.read_Mol2_file("DRCN5T-GSOpt-2.mol2")[0] ) # 432

M = [pc71, ]*266 + [drcn, ]*432
assert sum(map(lambda g:g.Atom_Count, M)) <= 100000, "Gro file does not support a system that contains more than 100,000 atoms"
random.shuffle(M)

B = Box((150, 150, 150))
RETRY = 1<<100

for ind, m in enumerate(M):

    while RETRY > 0:

        X, Y, Z = numpy.random.random(size=3)*B.size
        Rx, Ry, Rz = numpy.random.random(size=3)*2*numpy.array((numpy.pi, numpy.pi, numpy.pi))

        m.Rota = numpy.array([Rx, Ry, Rz])
        m.Tran = numpy.array([X, Y, Z])

        if B.Add_mol(m):
            print "INSERT Succ! - No. %d"%ind
            break

        else:
            RETRY -= 1
            print "RETRYING"


# A, B, Res, Box_size = B.Dump_box()

with open("t.mol2", "w") as f:
    print >>f, B.Box_2_mol_String("DRCN5T_PC71BM")
