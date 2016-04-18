# Installation in Linux/Mac #

  1. Checkout BCSLIB from http://code.google.com/p/bcslib/, which provides C++ support for the mex functions of this toolbox;
  1. Set environment variable _BCSLIB\_HOME_ in your system. For example, in Linux, you can add the following statement to your _.bashrc_ file:
```
export BCSLIB_HOME=<the root directory path of the BCSLIB>
```
  1. If you have not already done so, checkout _smitoolbox_
  1. Open MATLAB
  1. Go into the _tbman_ directory under the toolbox directory.
  1. Run _smi\_addpath_ to add paths to MATLAB. If you want to keep the path setting for future sessions, you may use the command _savepath_.
  1. Run _smi\_build\_mex_ to build all the mex files.
  1. Done.