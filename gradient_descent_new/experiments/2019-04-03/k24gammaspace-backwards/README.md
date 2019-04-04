Trying to generate gamma k24 phase plots for omega = 10 and Lambda = 1.0,10.0,100.0. In this folder, I'm scanning through gamma space backwards, to try and stay in the frustrated fibril phase as long as possible. Then comparing the scan forward and scan backward calculations, I can hopefully find the phase transition line and the critical point. To calculate phase plot (on compute canada), simply run:

$ sbatch run_cc.h 1.0 10.

$ sbatch run_cc.h 10.0 10.

$ sbatch run_cc.h 100.0 10.


for Lambda = 1.0, 10.0, or 100.0, respectively.