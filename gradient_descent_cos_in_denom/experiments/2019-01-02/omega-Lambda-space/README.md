THIS IS A FAILED ATTEMPT! SEE THE 2019-01-04 DIRECTORY FOR A NEWER VERSION WITH BETTER CODE.


This experiment calculated omega vs Lambda plots of the observables. For this run, I am looking at the space from omega = 0 to 30 and Lambda = 1 to 1000, in a 2d (3x3) grid of (gamma,k24) coordinates, with vertices (0.04,0.1), (0.04,0.9), (0.12,0.1), and (0.12,0.9).

To run this experiment, I need to submit several batch jobs to compute canada, by run "sh run_cc.h gamma k24 omega". Note that I'm NOT running it as sbatch, but allowing the bash script to submit it to the scheduler instead (for dependency reasons). runfirst.h and runsecond.h should not be run manually, only through run_cc.h should they be called.

The Lambda values to be calculated (i.e. Lambda = 0 to Lambda = 1000 in this case) must be set manually within the python script "scan2d.py". All other parameters (beside gamma, k24, omega, and Lambda range) are set in the input.dat file.