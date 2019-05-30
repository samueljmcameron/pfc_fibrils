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
    
    def __init__(self,readparams,tmp_path="../../../tmp_data/",
                 params=None,executable="../../../bin/double_well_scan"):

        self.readparams = readparams
        self.tmp_path = tmp_path
        self.executable = executable

        if (params == None):
            self.params = self.readparams.params
        else:
            self.params = params

        return



    def run_exe(self,valgrind = False):
        # run c executable to determine psi(r), R, delta, etc.

        if valgrind:
            subprocess.run(['valgrind','--track-origins=yes',self.executable,self.tmp_path,
                            *self.params.values()],check=True)
        else:
            subprocess.run([self.executable,self.tmp_path,*self.params.values()],check=True)

        return


    def write_suffix(self,p_list):
        
        suffix = '_'.join([f'{float(s):.4e}' for s in p_list])
        
        return suffix

    def mv_file(self,mname,newname=None):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.readparams.write_suffix()


        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        if newname != None:
            newfname = f"_{newname}_{suffix}.txt"
        else:
            newfname = fname

        mvto = f"data/{newfname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

