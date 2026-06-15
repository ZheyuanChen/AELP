# AELP (Arbitrary EPOCH Laser Profile)

## Obtaining a copy of modified EPOCH
In order to use the feature, you have to obtain my modified copy of EPOCH. You can do this by visiting [my fork on Github](https://github.com/ZheyuanChen/epoch_dev/tree/my-epoch-mods) or email me zheyuan.chen@york.ac.uk. If you choose the former, please download the branch "my-epoch-mods". This is VERY important as I don't think I keep the modification in the main branch. If you choose the latter, I will redirect you the the former, so you will use Github anyway.

Upon obtaining a copy, you can compile it in the usual way. 

## Feature introduction
I have modified the source code such that you are able to initialise an arbitrary laser profile in **2D simulations**. The work in 3D is under planning and the work in 1D is never intended. 

I think it's fairly simple to use. 

### Step 1: Generating a laser profile
You need to generate a file containing the desired laser profile. The file must be named "temporal_spatial_profile.dat" and must be contained in the same directory as the input.deck. I use the USE_DATA_DIRECTORY method to run EPOCH, and the modification is written assuming one uses this method. I haven't tested the compatibility with other methods. 

The tutorial contains a sample Python script on generating such a file (tutorials/arbitrary_laser_profile/generate_spatio_temporal_profile.py). 

The format of the file: 
First line: $N_t$ $N_y$ # the number of points along t and y. Note that I haven't enabled EPOCH to accept an injection on x-boundaries nor multiple laser beams.
Second line: y-coordinates
Third line: t-coordinates
Subsequent: the $N_t \cross N_y$ matrix of laser amplitude, normalised to 1. Each row represents a time slice, containing all y-coordinates.

Warning: I actually forgot how I wrote the source code, so the matrix may need to be transposed. I really need to check this.

There is a way to inject a spatial-only profile (and you can specify time profile using an analytical expression as before) but I have forgotten how to do this. 

## Step 2: Modify your input.deck
In the laser block, instead of initialising the profile by passing an analytical expression, write "profile = custom". There is a way to pass a spatial-only profile, but I forgot how:(

I have included some tutorial in this repo. You can have a look at them, particularly how the input.deck is written, and there are example scripts that generate spatial(-temporal) laser profiles for testing purposes.
