import numpy as np
import subprocess
import sys
import os

class MultiRun(object):
    # similar to singlerun.py, perform a single run of calculation for psi(r),
    # R, eta, ..., etc., except choose which parameters to vary by specifying
    # it when creating the object.

    # attributes:
    #  datfile - the data file that you read parameters from
    #  scan - a dictionary where you can pre-specify certain parameters (vs
    #         reading them in through the datfile).
    #  params  - an array that will have the string of parameters in it
    #  tmp_path - the path of where the output files are stored initially
    #  executable - the executable that creates and writes the output files
    
    def __init__(self,datfile,scan,params=None,tmp_path="../../../tmp_data/",
                 executable="../../../bin/gamma_k24_singlepoint"):

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

    def run_exe(self):
        # run c executable to determine psi(r), R, delta, etc.

        subprocess.run([self.executable,self.tmp_path,*self.params],check=True)

        return


    def write_suffix(self,p_list):
        
        suffix = '_'.join([f'{float(s):.4e}' for s in p_list])
        
        return suffix

    def get_xvals(self):
        
        suffix = self.write_suffix(self.params[:6])
        
        fname = f"_observables_{suffix}.txt"

        fullpath = f"data/{fname}"

        with open(fullpath,"r") as f:
            line = f.readline()
            E,R,eta,delta,surftwist = line.split()
        
        return R,eta,delta
    
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


    def concatenate_observables_with_direction(self,scan_dir):

        xval_params = [self.params[0],self.params[1],self.params[3],
                       self.params[5]]

        newfname = f"data/_observables_scan{scan_dir}_{self.write_suffix(xval_params)}.txt"

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
