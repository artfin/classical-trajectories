#!/usr/bin/python
import subprocess

programs = {
    'mpi-main': ('mpi-main', 3),
}

sys_call = '{0} -n {1} ./{2}'.format('mpirun', programs['mpi-main'][1], programs['mpi-main'][0])

subprocess.call([sys_call], shell = True)

