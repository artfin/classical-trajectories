#!/usr/bin/python
import subprocess

programs = {
    'main': ('main', 2),
}

sys_call = '{0} -n {1} ./{2}'.format('mpirun', programs['main'][1], programs['main'][0])

subprocess.call([sys_call], shell = True)

