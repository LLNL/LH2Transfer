
# LH2 Transfer Simulation

The code LH2TS (for “LH2 Transfer Simulation”) enables to calculate the wasted H2 due to transfer and boil-off while using LH2 as a fuel for transportation. More specifically, a reduced lumped-parameter model simulates thermodynamic states in 2 vessels (one feeding, one receiving) using mass and energy balances for each control volume (vapor and liquid, in each vessel) and relying on 2-phase behavior and heat transfer models specific to cryogenic fluids (self consistent theory of dynamical condensation/evaporation – see references for details).
This code was developed and tested using MATLAB R2013b and REFPROP DLL version 9.1 for parahydrogen [1], on a laptop PC running on 64-bit Windows 10 Enterprise (2016), with an Intel Core -7-6600U CPU @ 2.60 GHz and a 16.0 GB RAM.

# Description of the files

The simulation code consists of 9 MATLAB files that should be located in the same folder:
-	runNominal.m : runs the overall code

-	inputs_TrailerToDewar.m : initializes the inputs and boundary conditions

-	LH2Control.m : controls the different steps of the fill process

-	LH2Simulate.m : solves the governing ordinary integro-differential equations using an ODE solver (odes15s)

-	Data_extraction.m : post-processes output data

-	cylVtoH.m: function used in the post-process step. Converts volume of LH2 in a horizontal cylinder to a height of fluid

-	gasFlow.m: function used in the post-process step. Calculates the flow rates (chocked or not) between 2 vessels.

-	vaporpressure.m: function used in the post-process step. Calculates vapor pressure (2 phases or not) based on temperature and density.

-	plotLH2Data.m : plots output data, in multiple separate figures


Running runNominal.m with the MATLAB software will solve governing equations and save the results in “output.txt”, located in the same folder as the other files. Please refer to runNominal.m for the nomenclature of the output file.

The code architecture and underlying physical phenomena are based on the work of Muratov, Daigle, Osipov et al [2-5], who kindly shared their original code that was released as open-source. More information on the code can be found in the references. 

Modifications were made to the original version, including new shape factors, real gas equations of state linked to REFPROP for both thermodynamics and transport properties, new energy equations for the liquid phases, energy balances based on internal energies (as opposed to temperature). Polynomial expressions derived from REFPROP were used in some instances (e.g. saturated liquid density as a function of internal energy, enthalpy of vaporization as a function of film temperature).

The code has been verified only for “typical” thermodynamic states of H2 (i.e. temperatures, pressures and densities close to or within the 2 phase region). Atypical scenarios such as the cool-down with LH2 of a warm vessel have not been verified and the code would likely need additional modifications in order to adequately simulate those scenarios.

Please not that package only works if a REFPROP subroutine compatible with MATLAB is installed on the computer. For linking REFPROP with MATLAB, make sure to follow instructions that can be found online, such as http://trc.nist.gov/refprop/LINKING/Linking.htm (last accessed 12/7/2017). If running into issues with .dll from MATLAB, make sure you install REFPROP in “Program Files (x86)” and that REFPRP64_thunk_pcwin64.dll is in that directory.

# References

[1] Lemmon, E.W., Huber, M.L., McLinden, M.O.  NIST Standard Reference Database 23:  Reference Fluid Thermodynamic and Transport Properties-REFPROP, Version 9.1, National Institute of Standards and Technology, Standard Reference Data Program, Gaithersburg, 2013.

[2] V. V. Osipov and C. B. Muratov, “Dynamic condensation blocking in cryogenic refueling,” Applied Physics Letters, vol. 93, no. 22, pp. 224 105–1–4, 2008.

[3] V. Osipov, C. Muratov, M. Daigle, M. Foygel, V. Smelyanskiy, and A. Patterson-Hine, “A dynamical
physics model of nominal and faulty operational modes of propellant loading (liquid hydrogen): From space shuttle to future missions,” NASA Ames Research Center, Tech. Rep. NASA/TM-2010-216394, Jul. 2010.

[4] M. Daigle, S. Lee, C. Muratov, S. Osipov, V. Smelyanskiy, and  Michael Foygel, “LH2 Matlab Simulation User Manual”, December 3 2010

[5] M. Daigle, M. Foygel and V. Smelyanskiy, “Model-based diagnostics for propellant loading systems”,  2011 IEEE Aerospace Conference, 5-12 March 2011

# License

LH2TS is released under a NASA open source agreement license, version 1.3. For more details see the [LICENSE](https://github.com/LLNL/LH2Transfer/blob/master/LICENSE.md) file

LLNL-CODE-748554
