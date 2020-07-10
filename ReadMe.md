![](https://raw.githubusercontent.com/Novarizark/cesm-cam-forcing-mod/master/forcing_module.jpg)

**CESM CAM FORCING MODULE (CCFM)** is an open source modification module for the [CESM](http://www.cesm.ucar.edu/) users using external forcing files to conduct sensitive experiments in the CAM workflow. This flexible architecture enables you to deploy the module to a series of different versions of CESM (with careful check), although the source code was developed based on CESM1.2.2.

CCFM was originally developed by LZN to conduct prescribed convective heating forcing experiment in 2015. Then, he found that the similar technique was general enough to be applicable in a wide variety of other questions by modifying the model prognostic physical tendencies. Some complex jobs e.g. Nudging, can also be implemented with the similar basis. Therefore, he decided to develop a common module embedded into the CESM architecture, which can be simply used by other scholars to conduct related experiments.

## Installation

CCFM is embedded into the CESM source code, to install the current release of CCFM, you need first download this repo to your targeted machine.

If you have installed `git` on a Linux platform, simply type:
```
git clone git@github.com:Novarizark/cesm-cam-forcing-mod.git 
```

Alternatively use `wget`:
```
wget https://github.com/Novarizark/cesm-cam-forcing-mod/archive/master.zip
unzip master.zip
```

Next, copy the modified `namelist_definition.xml` to cover the original one in CAM:
```
cd cesm-cam-forcing-mod
cp namelist_definition.xml $CCSM_ROOT/models/atm/cam/bld/namelist_files/
```
`$CCSM_ROOT` is the root directory of the CESM. 

Finally, after setup the case, copy the `src.cam` to the case SourceMods directory.
```
cp -r src.cam $YOUR_CASE_ROOT/SourceMods/
```

## Usage

* `cam_forcing_mod.F90` is the key module that define, register, initial and execute the forcing into the CAM workflow.

Last Updated: April 30, 2019
