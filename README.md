Dark matter capture on electrons in White Dwarfs

================================================================
Code Compilation
================================================================

Use:
First complie the code with make, then run the run_capelec.exe with the EoS number as input i.e.

    make
    ./run_capelec

which runs the code for the default 1.38367 M_sun WD on half the system cores.

OR 

just run 
    ./comp_cap
for the same resuls.

There is also a python script to do the same thing, but can take the input 'log' to log terminal outputs to /logs/term_err.txt and /logs/term_out.txt

    python3 py_run.py log

================================================================
Using the C code
================================================================

Use the routines in the main file to generate the relevent profiles:

    for all operators:
        input: Log10(m_chi_MIN), Log10(m_chi_MAX), step size

    for single operators:
        input: Operator number (11 for constant XS), Log10(m_chi_MIN), Log10(m_chi_MAX), step size

If running on Spartan (or other cluster with OpenMP), make sure the variable USE_SPARTAN is set to 1.
Else set to 0, with default number of cores set to half the system maximum. 

On spartan, first submit the "remake_exec.slurm" to make sure the right executable is there. Then submit the "job_WD.slurm" which creates seperate jobs for each  EoS configuration. This can be changed in the job file.

Multiscattering works up to 10^9 MeV. Need to change the integrator to the approximated 1D integral, that should increase the range. 