function f = pewarp_fminunc_regularised_costwrapper(x)

global bX bY dbY bZ target8 source mask krn xscale regstrength

f = mex_pewarpcost_regularised(single(x), bX, bY, dbY, bZ, target8, source, mask, krn, xscale, regstrength);
