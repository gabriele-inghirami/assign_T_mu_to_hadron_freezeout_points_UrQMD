## Quick infomation

This repository contains stuff used for: https://arxiv.org/abs/2007.06440.
The repository has been assembled two years after the study, just collecting code left here and there, without description of the workflow and with poorly documented scripts.
This is not the right way to proceed, in more recent works code and workflow have been archived in a more clear format immediately after the submission of the article,
but it is still better than throw away the files and leave nothing...


Main script:

get_info_chemfreeze_v2.2.1.py - it associates the temperature, the baryon chemical potential and the transverse velocity obtained with make_cg_plots_v2.2.0.py to a list of hadron coordinates. It is possible to restore the full format, which includes also the distributions with respect to transverse momentum and the rapidity of the hadrons. The script creates an extended list of hadrons in ascii format and several kind of histograms.

In assorted_plotting_scripts:

- combine_info_v2.2.0.py - it combines the histograms obtained with get_info_chemfreeze_v2.2.1.py from different (but homogeneous) datasets into a single file

- make_plots_v2.4.0.py - it plots the histograms computed by get_info_chemfreeze_v2.2.1.py (for a single collision energy)

- make_time_plots.py - it prepares time histograms (but it does not make the plots) from the output of get_info_chemfreeze_v2.2.1.py

- make_dNdT_dNdmu_plots_noa.py, make_dNdT_dNdmu_plots.py - it plots these distributions, prepared by the get_info_chemfreeze_v2.2.1.py, to compare simultaneously the results for different collision energies

- make_dvar_dtime_plots_noa.py, make_dvar_dtime_plots.py - it plots these distributions, prepared by the make_time_plots.py, to compare simultaneously the results for different collision energies

- quick_avg_min_temp.py - it computes some averages from the ascii files created by get_info_chemfreeze_v2.2.1.py
