# Nordic5
* Contact: Joakim Björk, joakbj@kth.se

This repository presents a 5 machine model of the Nordic synchronous grid, developed in Matlab Simulink Simscape Electrical. For full details of the model and its design choices, see the document **[Nordic5_doc](Nordic5_doc.pdf)**.

# Dependencies
1. Matlab (R2020b)
2. Simulink 
3. Simscape Electrical

# How to use model

Most of the implemented functions use the library file **Nordic32_Lib.slx**, wind turbine data in the **Wind_NREL_data**-folder and functions in the **Functions**-folder. Run **open_path_here.m** to load all dependencies. Projects are sorted into different folders (see below). For all of these, running **main_.m** allows the user to recreate results, either from simulation, or from saved examples.

## Base_Case_Nordic5_5000MW

- Main File: **main_PSS_controller_design.m**

These files are used to desing different load flows and to tune power system stabilizers. Current version is documented in **Nordic5_doc.pdf**

## Ensemble_Nordic5_5000MW, Ensemble_two_hydro_DVPP

Contains the simulation examples used in [[1]](#1).

Main Files: 

- **main_Ensemble_Nordic5_ideal.m** 
- **main_Ensemble_Nordic5_wind_hydro.m**
- **main_ensemble_two_hydro_DVPP.m**

## Wind_and_hydro_DVPP, Wind_Nordic5_5000MW, Wind_open_loop_step

Contains the simulation examples used in [[2]](#2). The implemented wind turbine model is a modified version of the NREL 5 MW baseline wind turbine model [[3]](#3).

Main Files: 

- **main_wind_and_hydro_DVPP.m** 
- **main_Wind_Nordic_5000MW_110.m**
- **main_wind_open_loop.m**

## References
<a id="1">[1]</a> 
J. Björk, K. H. Johansson, and F. Dörfler, “Dynamic virtual power plant design for fast frequency reserves: Coordinating hydro and wind,” IEEE Control Netw. Syst., 2022. 
<br />
https://doi.org/10.1109/TCNS.2022.3181553


<a id="2">[2]</a> 
J. Björk, D. V. Pombo, and K. H. Johansson, “Variable-speed wind turbine control designed for coordinated fast frequency reserves,” IEEE Trans. Power Syst., vol. 37, no. 2, pp. 1471 – 1481, Mar. 2022. 
<br />
https://doi.org/10.1109/TPWRS.2021.3104905
<br />
https://arxiv.org/abs/2108.02427


<a id="3">[3]</a> 
J. Jonkman, S. Butterfield, W. Musial, and G. Scott, “Definition of a 5-MW reference wind turbine for offshore system development,” NREL, USA, Tech. Rep., 2009.
