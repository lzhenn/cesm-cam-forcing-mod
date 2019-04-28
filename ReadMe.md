![](https://raw.githubusercontent.com/Novarizark/cesm-cam-forcing-mod/master/forcing_module.jpg)

**CESM CAM FORCING MODULE (CCFM)** is an open source modification module for the [CESM](http://www.cesm.ucar.edu/) users using external forcing files to conduct sensitive experiments in the CAM workflow. This flexible architecture enables you to deploy the module to a series of different versions of CESM, although the source code was developed based on CESM1.2.2.

CCFM was originally developed by LZN to conduct prescribed convective heating forcing experiment in 2015. Then, he found that the similar technique was general enough to be applicable in a wide variety of other questions by modifying the model prognostic physical tendencies. Some complex jobs e.g. Nudging, can also be implemented with the similar basis. Therefore, he decided to develop a common module embedded into the CESM architecture, which can be simply used by other scholars to conduct related experiments.

Last Updated: April 29, 2019
