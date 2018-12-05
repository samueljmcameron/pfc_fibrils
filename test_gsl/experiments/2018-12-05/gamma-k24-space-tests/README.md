This experiment aims replicate a previous experiment (from 2018-11-28/gamma-k24-at0.15-0.5) to show how turning omega and Lambda on affects the phase transition behaviour of R, eta, and delta at a specific value of (gamma,k24) = (0.15,0.5).

I just changed the convergence criterion from |dEdy|<CONV_MIN to |dEdy|<CONV_MIN*(1+|E|) in the code and wanted to make sure nothing changed.