A few pointers:

1) After running aamod_coregister_unwarp the functional and structural images are coregistered, so you can directly feed them into normalisation.

2) It is assumed that the structural has been brain extracted (with BET, BSE, etc.), and ideally bias corrected. The mean EPI is bias corrected within the script

3) This has only been tested with SPM8.

4) mex_pewarpcost_regularised is a compiled function that is built for Matlab 7.8. It will probably not work in older versions. If you want to recompile you can use the buildC99ompblas script, but it's probably quite a hassle as you'll need a recent version of gcc. (Summary: Just use matlab 7.8 or newer)

If you run into any problems, let us know!
Tool developer: Eelke (eelke.visser@donders.ru.nl)
Ported to AA: Alejandro (a.vicente.grab@gmail.com)

If you use this tool, please cite Eelke's work:
EPI DISTORTION CORRECTION BY CONSTRAINED NONLINEAR COREGISTRATION IMPROVES GROUP FMRI
E. Visser1,2, S. Qin1,3, and M. P. Zwiers1,2
1Donders Institute for Brain, Cognition and Behaviour, Radboud University Nijmegen, Nijmegen, Netherlands, 2Department of Psychiatry, Radboud
University Nijmegen Medical Centre, Nijmegen, Netherlands, 3Department of Neurology, Radboud University Nijmegen Medical Centre, Nijmegen,
Netherlands
Proc. Intl. Soc. Mag. Reson. Med. 18 (2010)