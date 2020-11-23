function pCanon = readlink(pIn)
p = java.io.File(pIn);
pCanon = char(p.getCanonicalPath());

