# Driving behavior assessment

A repository of an on course project. This project attempts to assess driving behavior impact on electric vehicles from physiological signals.

## Files for computing features from Physiological signals

"read.m" contains the lines for computign these feautres. ZCR and emd are functions called by the aforementioned.

## Files for image processing for feature extraction from vehicle data
 
 This code is conditioned to videos with 1280x720 resolution.
 
"perfil.m" contains the lines for driving profile, braking and acceleration events extraction depending on the video inserted.
"brake.m" contains the code for extracting the accelerate and brake pressure for an inserted video.

.mat files are necessary to be into the same location of "perfil.m" to make this code functional.
