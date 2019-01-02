# accept: serial, HOMO(in ev), LUMO(in ev)
import numpy, sys

D = numpy.genfromtxt(sys.argv[1])[:, 1:]
HOMOs = D[:, 0]
LUMOs = D[:, 1]
Ec = numpy.min(LUMOs)
Ev = numpy.max(HOMOs)

Nc = 2*sum( numpy.exp( (HOMOs - Ec)/-0.026 ) )
Nv = 2*sum( numpy.exp( (Ev - LUMOs)/-0.026 ) )

print "Nc = %3.2e"%Nc
print "Nv = %3.2e"%Nv