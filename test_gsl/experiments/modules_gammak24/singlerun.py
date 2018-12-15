# Last updated on December 5th, 2018. I merged the original multirun.py module with
# the singlerun.py module, as it was unnecessary to have two modules which essentially
# do the same thing.


import numpy as np
import subprocess
import sys
import os

class SingleRun(object):
    # Perform a single run of calculation for psi(r), R, eta, ..., etc.
    # You can choose which parameters to vary by specifying a dictionary
    # of values "scan" when creating the object.

    # attributes:
    #  datfile - the data file that you read parameters from
    #  scan - a dictionary where you can pre-specify certain parameters (vs
    #         reading them in through the datfile).
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,readparams,tmp_path="../../../tmp_data/",
                 params=None,executable="../../../bin/gamma_k24_singlepoint"):

        self.readparams = readparams
        self.tmp_path = tmp_path
        self.executable = executable
        if (params == None):
            self.params = self.readparams.params
        else:
            self.params = params

        return

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        subprocess.run([self.executable,self.tmp_path,*self.params.values()],check=True)

        return

    def get_xvals(self):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.readparams.write_suffix()
        
        fname = f"_observables_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            E,R,eta,delta,surftwist = line.split()
        
        return R,eta,delta
    
    def mv_file(self,mname):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.readparams.write_suffix()

        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        mvto = f"data/{fname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

    def add_datastring(self,vars):
        
        a = [float(self.params[var]) for var in vars]

        return '\t'.join(map("{:13.6e}".format,a))

    def concatenate_observables(self,vars,scan_dir=''):

        suffix1 = self.readparams.write_suffix(suffix_type="save")

        newfname = f"data/_observables_{scan_dir}_{suffix1}.txt"

        with open(newfname,"a+") as f1:
    
            suffix2 = self.readparams.write_suffix()

            fname = f"data/_observables_{suffix2}.txt"

            with open(fname) as f2:

                for line in f2:

                    f1.write(f"{self.add_datastring(vars)}\t{line}")

            if os.path.isfile(fname):

                os.remove(fname)
        
        return
