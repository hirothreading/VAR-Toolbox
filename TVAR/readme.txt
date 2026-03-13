The code runs in this sequence:

(1) estimatenew1.m. This estimates the reduced form model. Note that options for the code are set up to line 43 where the data is loaded.

(2)getirf1SD_r1.m: IRF in regime 1. Currently uses cholesky decomposition (see line 36 and 37)

(3)getirf1SD_r2.m: same as above but for regime 2


You can plot the results for this example by then running plotsupply_all_regime


This model is described in https://ideas.repec.org/a/red/issued/14-103.html