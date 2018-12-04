import numpy as np
import subprocess
import sys
import os

class SingleRun(object):
    # perform a single run of calculation for psi(r), R, eta, ...

    # attributes:
    #  datfile - the data file that you read parameters from
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,datfile,params=None,tmp_path="../../../tmp_data/",
                 executable="../../../bin/gamma_k24_singlepoint"):

        self.datfile = datfile
        self.tmp_path = tmp_path
        self.executable = executable

        if (params == None):
            self.params = self.read_params()
        else:
            self.params = params
        
        return

    def read_params(self):
        # read input parameters from datfile.

        params = []
        
        with open(self.datfile) as f:

            for line in f:

                (key,val) = line.split()
                params.append(val)

        return params

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        subprocess.run([self.executable,self.tmp_path,*self.params],check=True)

        return


    def write_suffix(self,p_list):
        
        suffix = '_'.join([f'{float(s):.4e}' for s in p_list])
        
        return suffix

    def mv_file(self,mname):
        # move a file from the temporary file folder to the data folder.
    
        suffix = self.write_suffix(self.params[:6])

        fname = f"_{mname}_{suffix}.txt"

        mvfrom = f"{self.tmp_path}{fname}"

        mvto = f"data/{fname}"

        subprocess.run(["mv",mvfrom,mvto],check=True)

        return

    def concatenate_observables(self):

        xval_params = [self.params[0],self.params[1],self.params[3],
                       self.params[5]]

        newfname = f"data/_observables_{self.write_suffix(xval_params)}.txt"

        with open(newfname,"a+") as f1:
    
            fname = f"data/_observables_{self.write_suffix(self.params[:6])}.txt"

            with open(fname) as f2:
            
                for line in f2:
                    add_omega_Lambda = '\t'.join([f'{float(self.params[2]):13.6e}',
                                                  f'{float(self.params[4]):13.6e}'])
                    f1.write(f"{add_omega_Lambda}\t{line}")
            if os.path.isfile(fname):
                os.remove(fname)
            
        return
