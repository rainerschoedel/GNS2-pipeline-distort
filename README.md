# GNS2-pipeline-distort
GNS2 pipeline based on pre-reduced and distortion-corrected data provided by Hervé Bouy 

# Initial commit on 26 January 2025.

Execute the code in the following order

The pre-reduced and distortion corrected images provided by Hervé must be located in "../cubes/" and the
corresponding weights in "../weights/".


# (1) bgcorrect.py

Find and correct sky offset (there may be negativities because of Hervé's pipeline.
Application of this code is optional and typically not necessary.

# (2) lxp.py

Create a long exposure mosaic image.

# (3) findstars.pro

Find stars in long expoure image. Do this for the four chips to be able to run the code in parallel.

# (4) subcubes.py

Create sub-cubes.

# (5) sort_subcubes.py

Sort the subcubes so that they can be used correctly by holo_full.pro

# (6) vvv_stars.py

Find astrometric reference stars in VVV_X data.

# (5) align_VVV.py

Find preliminary alignment stars for each chip; to be used later in calibrate.pro (photometry).

# (6) runholo.pro

Run holography on the subcubes.

