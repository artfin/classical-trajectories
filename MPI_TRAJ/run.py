#!/usr/bin/python

import sys
import subprocess

instances = 4 

programs = {
    'main': ('main', instances),
    'mpi_reader': ('mpi_reader', instances),
}

program_to_run = sys.argv[1] if len(sys.argv) > 1 else None
if not program_to_run in programs:
    print 'Must enter programe to run. Possible programs are: {0}'.format(programs.keys())
else:
    sys_call = '{0} -n {1} ./{2}'.format(
        'mpirun', programs[program_to_run][1], programs[program_to_run][0])
    
    if len(sys.argv) > 2:
        for arg in sys.argv[2:]:
            sys_call += (" " + arg) 
    
    print(sys_call)

    subprocess.call( [sys_call], shell = True )




