import numpy
from mpi4py import MPI

comm = MPI.COMM_WORLD

print("hello", comm.rank)

comm.Barrier()

