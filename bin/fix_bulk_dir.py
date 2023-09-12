#!/usr/bin/python

import sys
import os

bulk_outdir = sys.argv[1]
bulk_name = sys.argv[2]
sc_name = sys.argv[3]

assert os.path.exists(bulk_outdir)

if os.path.exists(os.path.join(bulk_outdir, sc_name)):
    print(f'{sc_name} folder alredy present...')
elif not os.path.exists(os.path.join(bulk_outdir, sc_name)) and os.path.exists(os.path.join(bulk_outdir, bulk_name)):
    os.chdir(bulk_outdir)
    os.system(f'ln -s {bulk_name} ./{sc_name}')
    print(f'{sc_name} folder not present. Found bulk one that map to it.')
else:
    print(f'Neither {sc_name} nor {bulk_name} are present... Skip this one!')
