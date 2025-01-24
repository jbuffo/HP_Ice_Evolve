# HP_Ice_Evolve
1D Thermophysical Simulations of Layered Ice Sheets (Including High Pressure Ices and Their Thermodynamics)

THIS IS A PLACEHOLDER README - I WILL UPDATE WITH A MORE COMPLETE DESCRIPTION OF THE CODE SOON

HP Ice Evolve

This code was designed for and accompanies the manuscript “Liquid water on cold exo-Earths via basal melting of ice sheets” (https://www.nature.com/articles/s41467-022-35187-4) and simulates the One Dimensional thermophysical evolution of very thick ice sheets - that can reach thicknesses and/or pressures that induce the transition to High Pressure (HP) Ice phases.

This code makes use of another open-source thermodynamic code called SeaFreeze (https://github.com/Bjournaux/SeaFreeze). You will need to download the MATLAB functions of SeaFreeze to run HP Ice Evolve (this code).

The code is initiated by describing the dimensions, thermal, and pressure environments of an ice sheet. An initial pressure, adiabatic temperature, and ice phase profile will be calculated. The code will then evolve the ice sheet given top and bottom thermal boundary conditions, calculating phase changes, evolving the conductive/convective thermal profile, and simulating any meltwater generation and transport through or beneath the ice layers.

The basic user input variations will be carried out by modifying lines in the 1-100 region of the HP_Ice_Evolve_Convect_v2.m file. This is also the file that will be run to execute the simulation.

Please feel free to reach out to me with any questions/comments/etc. at jacob.j.buffo@dartmouth.edu and I will try to get back to you as soon as possible - however, apologies in advance for any slow replies as work keeps me busy and I am not currently funded to maintain/update this code.

Best,

jacob
