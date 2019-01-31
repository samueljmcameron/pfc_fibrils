import numpy as np
import subprocess
import sys
import os

class DoubleWell(object):
    # perform a single run of calculation for psi(r), R, eta, ...

    # attributes:
    #  datfile - the data file that you read parameters from
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,datfile='data/input.dat',scan={},params=None,
                 tmp_path="../../../tmp_data/",
                 executable="../../../bin/double_well_scan"):

        self.datfile = datfile
        self.tmp_path = tmp_path
        self.executable = executable
        self.scan = scan
        if (params == None):
            self.params = self.read_params()
        else:
            self.params = params
        
        return

    def set_param(self,key,val):
        
        for scan_key,scan_val in self.scan.items():

            if (key == scan_key):
                return scan_val

        return val

    
    def read_params(self):
        # read input parameters from datfile.

        params = []
        
        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()
                params.append(self.set_param(key,val))

        return params


    def run_exe(self,valgrind = False):
        # run c executable to determine psi(r), R, delta, etc.

        if valgrind:
            subprocess.run(['valgrind','--track-origins=yes',self.executable,self.tmp_path,
                            *self.params],check=True)
        else:
            subprocess.run([self.executable,self.tmp_path,*self.params],check=True)

        return


    def write_suffix(self,p_list):
        
        suffix = '_'.join([f'{float(s):.4e}' for s in p_list])
        
        return suffix

    def mv_file(self):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.write_suffix(self.params[:6])

        fname = f"_Evst_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        mvto = f"data/{fname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return
