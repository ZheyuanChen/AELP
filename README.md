# AELP (Arbitrary EPOCH Laser Profile)

## Update (30/06)
File format overhaul -- stick to EPOCH's convention (binary files)
Several bug fixing
Reformating epoch-mod to conform to EPOCH's coding style requirements
Warning: due to the file format change, every input.deck and .dat files written before 30/06 won't work

## The Very Latest Update (25/06)
Phase injection implemented and tested in 2D -- ready to use. 
LASY integration in progress.
3D implementation postponed.

## THE Latest Update (24/06)
EPOCH3d spatial profile ready to use. Spatial-temporal profile underconstruction. 
LASY integration under construction.
A potential overhaul on the file loader (changing from dat files to binary files or keeping both formats).

## Latest Update (19/06)
Laser profile injection in EPOCH2d is stable and ready to use. Tests have been done and show good results.
EPOCH3d modification is under construction.

## Obtaining a copy of modified EPOCH
In order to use the feature, you have to obtain my modified copy of EPOCH. You can do this by visiting [my fork on Github](https://github.com/ZheyuanChen/epoch_dev/tree/my-epoch-mods) or email me zheyuan.chen@york.ac.uk. If you choose the former, please download the branch "my-epoch-mods". This is VERY important as I don't think I keep the modification in the main branch. If you choose the latter, I will redirect you the the former, so you will use Github anyway.

Upon obtaining a copy, you can compile it in the usual way. 

## Installation
If you already have a working virtual environment, you can go straight to Jupyter notebooks and play them yourself. If you wish to create a new venv for this project, you can download all required Python package by `uv pip install -e .` assuming you are at the top directory of this project and you use uv. 

## Feature introduction (OUTDATED -- Documentation Pending update)
I have modified the source code such that you are able to initialise an arbitrary laser profile in **2D simulations**. The work in 3D is under planning and the work in 1D is never intended. 

I think it's fairly simple to use. 

### Step 1: Generating a laser profile
You have three options to initialise a laser profile.
1. The usual way: in this case the only thing you need to do is add a line "use_custom_profile = F" in your laser block and proceed as usual.
2. Use a customised file specifiying the spatial profile only: you need to generate a file named "spatial_profile.dat". In the first line you need to write the total number of points (e.g. 500). In the following lines, you write the spatial coordinate (in metres) and the laser amplitude (normalised to 1), separate by a space (e.g. -10.0e-6 0.5). Note that you can still pass an analytical function as the temporal profile using the usual way (t_profile = gauss(.,.,.) for example)
3. Use an array specifying the temporal-spatial profile: yuo need to generate a file named "temporal_spatial_profile.dat". See below for its format as it's a little bit complicated.

The format of the file "temporal_spatial_profile.dat": 
First line: $N_t$ $N_y$ # the number of points along t and y. Note that I haven't enabled EPOCH to accept an injection on x-boundaries nor multiple laser beams.
Second line: y-coordinates
Third line: t-coordinates
Subsequent: the $N_t \cross N_y$ matrix of laser amplitude, normalised to 1. Each row represents a time slice, containing all y-coordinates.

Warning: I actually forgot how I wrote the source code, so the matrix may need to be transposed. I really need to check this.

### Step 2: Modifying your input.deck
There are two newly-added logical flags that you need to specify in the laser block in input.deck. They are "use_custom_profile" and "use_spatiotemporal_profile". 
If you wish to use a custom laser profile, set "use_custom_profile = T". The default is F. If the flag is set to F, you can set the profile in the usual way (by passing an analytical expression).

On top of this, if you simply want to set a discrete spatial profile (i.e. E = E(y)), set "use_spatiotemporal_profile = F". The default is T. In this case, you need to include a file named "spatial_profile.dat" in the same directory as your input.deck. Note that I use and have only tested the "USE_DATA_DIRECTORY" method for running EPOCH (i.e. create a file named USE_DATA_DIRECTORY in epoch/epoch2d containing the input.deck's directory and call EPOCH without passing any directory). 

If you want to use a discrete temporal-spatio profile (i.e. E = E(t,y)), set "use_spatiotemporal_profile = T" and include a file named "temporal_spatial_profile.dat" in the same directory as your input.deck. 

## Known Issues:
Moved to Issues

## A Better Abbreviation:
AEPI(I) (Arbitrary EPOCH Profile Initialisation? Interface?) pronounced yippee

CELP (Customised/Customisation of/Customising EPOCH Laser Profile) pronounced in a similar way as Celtics


## Warnings
Not at the moment

## Developer's Section
TODO:
1. Maybe name things in a better and more consistent way.
2. epoch3d spatial-temporal profile
3. LASY integration
