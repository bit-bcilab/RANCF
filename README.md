# RANCF

The codes are implemented based on the [GOT-10k Matlab toolkit](https://github.com/got-10k/toolkit-matlab). 

### Requirments
* Matlab
* Nvidia GPU
* MatConvNet (need Visual Studio 2015 for compiling)

## instructions

#### Install and compile MatConvNet
Please refer to https://www.vlfeat.org/matconvnet/install/ for
installing and compiling. And Don't forget add its path in Matlab!

#### Test
Similarly, the path should be added in Matlab.

The tracker is implemented in the *got10k/trackers/* dir.

Other utility codes are in the *RANCF/* dir.

Finally, modify the settings in *quick_examples.m* and run it.
