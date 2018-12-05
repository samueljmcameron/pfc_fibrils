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
    
    def __init__(self,datfile,scan={},params=None,tmp_path="../../../tmp_data/",
                 executable="../../../bin/gamma_k24_singlepoint",
                 suffixlist=["K_{33}","k_{24}","\\Lambda","d_0",
                             "\\omega","\\gamma_s"]):

        self.datfile = datfile
        self.tmp_path = tmp_path
        self.executable = executable
        self.scan = scan
        self.suffixlist = suffixlist
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
        # from there.
        
        params = []
        
        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()

                params.append(self.set_param(key,val))

        print(params)

        return params

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        subprocess.run([self.executable,self.tmp_path,*self.params],check=True)

        return

    def make_plist(self):

        plist = []

        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()
                for p in self.suffixlist:
                    if p == key:
                        plist.append(self.set_param(key,val))

        print(plist)
        return plist


    def write_suffix(self,plist=None):
        
        if plist == None:
            plist = self.make_plist()

        suffix = '_'.join([f'{float(s):.4e}' for s in plist])
        
        return suffix

    def get_xvals(self):
        
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

    def concatenate_observables(self,scan_dir=''):

        xval_params = [self.params[0],self.params[1],self.params[3],
                       self.params[5]]

        newfname = f"data/_observables_{scan_dir}_{self.write_suffix(xval_params)}.txt"

        with open(newfname,"a+") as f1:
    
            fname = f"data/_observables_{self.write_suffix()}.txt"

            with open(fname) as f2:
            
                for line in f2:
                    add_omega_Lambda = '\t'.join([f'{float(self.params[2]):13.6e}',
                                                  f'{float(self.params[4]):13.6e}'])
                    f1.write(f"{add_omega_Lambda}\t{line}")
            if os.path.isfile(fname):
                os.remove(fname)
        
        return
