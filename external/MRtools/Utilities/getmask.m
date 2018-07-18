function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = true(size(y1));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end
