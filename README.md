# AELP (Arbitrary EPOCH Laser Profile)

## Obtaining a copy of modified EPOCH
In order to use the feature, you have to obtain my modified copy of EPOCH. You can do this by visiting [my fork on Github](https://github.com/ZheyuanChen/epoch_dev/tree/my-epoch-mods) or email me zheyuan.chen@york.ac.uk. If you choose the former, please download the branch "my-epoch-mods". This is VERY important as I don't think I keep the modification in the main branch. If you choose the latter, I will redirect you the the former, so you will use Github anyway.

Upon obtaining a copy, you can compile it in the usual way. 

## Feature introduction
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

### Step 2: Modify your input.deck
There are two newly-added logical flags that you need to specify in the laser block in input.deck. They are "use_custom_profile" and "use_spatiotemporal_profile". 
If you wish to use a custom laser profile, set "use_custom_profile = T". The default is F. If the flag is set to F, you can set the profile in the usual way (by passing an analytical expression).

On top of this, if you simply want to set a discrete spatial profile (i.e. E = E(y)), set "use_spatiotemporal_profile = F". The default is T. In this case, you need to include a file named "spatial_profile.dat" in the same directory as your input.deck. Note that I use and have only tested the "USE_DATA_DIRECTORY" method for running EPOCH (i.e. create a file named USE_DATA_DIRECTORY in epoch/epoch2d containing the input.deck's directory and call EPOCH without passing any directory). 

If you want to use a discrete temporal-spatio profile (i.e. E = E(t,y)), set "use_spatiotemporal_profile = T" and include a file named "temporal_spatial_profile.dat" in the same directory as your input.deck. 

## Known Issues:
1. If not using profile customisation, a weird bug is causing EPOCH to input zero fields if the analytical expression is time-independent. For example, profile = gauss(time,2*femto,2*femto) works, but profile = 1 doesn't. It doesn't matter if you write anything about t_profile. If there is any time dependency in profile, EPOCH will work properly and multiply the profile with t_profile. 

## Warnings
1. I only modify the code for lasers attached to the x_min boundary, so please at this stage only try to inject a laser from the left boundary. Also, I actually haven't tested what happens if you inject more than 1 lasers.
2. I use bilinear interpolation between the discrete points. Rigorous numerical testing is needed.

## Developer's Section
TODO:
1. Maybe name things in a better and more consistent way.
2. epoch3d
3. All boundaries
4. Testing
