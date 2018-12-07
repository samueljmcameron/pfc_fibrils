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
    
    def __init__(self,datfile="data/input.dat",scan={},tmp_path="../../../tmp_data/",
                 params=None,executable="../../../bin/gamma_k24_singlepoint",
                 loadsuf=["K_{33}","k_{24}","\\Lambda","d_0","\\omega","\\gamma_s"],
                 savesuf=["K_{33}","k_{24}","d_0","\\gamma_s"]):

        self.datfile = datfile
        self.tmp_path = tmp_path
        self.executable = executable
        self.scan = scan
        self.loadsuf = loadsuf
        self.savesuf = savesuf
        if (params == None):
            self.params = self.read_params()
        else:
            self.params = params

        return

    def set_param(self,key,val):
        # if the key passed as an argument here is already a key in the
        # self.scan dictionary, return the value in self.scan[key].
        # Otherwise, return val.
        
        for scan_key,scan_val in self.scan.items():

            if (key == scan_key):
                return scan_val

        return val
        

    def read_params(self):
        # read input parameters from datfile, unless the input parameter
        # is specified in self.scan dictionary, in which case, read it
        # from there. return a dictionary with parameter labels and names.
        # NOTE since this is running python 3.7 (or later), dicts are ordered.
        
        params = {}
        
        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()

                params[key]=self.set_param(key,val)

        return params

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        subprocess.run([self.executable,self.tmp_path,*self.params.values()],check=True)

        return

    def write_suffix(self,suffix_type="load"):
        # write the ending "suffix" of the file name (all of the trailing parameter values
        # in the filename. e.g. if fname = "energy_3.00000e+00_5.32000e-01.txt", then the
        # suffix="3.00000e+00_5.32000e-01", WITHOUT the ".txt".

        if (suffix_type=="save"):
            cpylist = self.savesuf
        else:
            cpylist = self.loadsuf

        plist = [float(self.params[s]) for s in cpylist]

        suffix = '_'.join([f'{float(s):.4e}' for s in plist])
        
        return suffix

    def get_xvals(self):
        # get the values in the x array, from the data file that HAS ALREADY BEEN MOVED FROM
        # tmp_data folder to true data folder, and BEFORE concatenation.

        suffix = self.write_suffix()
        
        fname = f"_observables_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            E,R,eta,delta,surftwist = line.split()
        
        return R,eta,delta
    
    def mv_file(self,mname):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.write_suffix()

        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        mvto = f"data/{fname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

    def add_datastring(self,vars):
        
        a = [float(self.params[var]) for var in vars]

        return '\t'.join(map("{:13.6e}".format,a))

    def concatenate_observables(self,scan_dir=''):

        suffix_type="save"

        newfname = f"data/_observables_{scan_dir}_{self.write_suffix(suffix_type=suffix_type)}.txt"


        # determine which variables are not shared between loadsuf and
        # savesuf, those will be the variables whose values are added
        # to the concatenated output file.
        vars = list(set(self.loadsuf).difference(self.savesuf))
        

        with open(newfname,"a+") as f1:
    
            fname = f"data/_observables_{self.write_suffix()}.txt"

            with open(fname) as f2:

                for line in f2:

                    f1.write(f"{self.add_datastring(vars)}\t{line}")

            if os.path.isfile(fname):

                os.remove(fname)
        
        return
