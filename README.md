# EyeRay

The presented code is an extension of the work developed by Gagnon et al. (App. Opt. 2014). In this method, a two-dimensional ray tracing procedure has been developed for an arbitrary number of surfaces and arbitrary surface shapes. The Liou and Brennan anatomically accurate eye model has been adapted and used for evaluating the method.

In order to run this code you must add Chebfun toolbox to your Matlab directory (http://www.chebfun.org/). Since Chebfun 5.2 does not provide arctangent for chebfun2 you can creat it manually or just add the file atan provided here to @chebfun2 folder.

Finally, run the code EyeRay.m to perform ray tracing over Liou and Brennan adapated model and EyeRayShow.m and/or gPSF.m to present the results graphically.

Comments are welcome
