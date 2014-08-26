function Wall=Wcons(W)

% assumes that W = W1111 W2222 W3333 W1112 W1113
%                  W1222 W2223 W1333 W2333 W1122
%                  W1133 W2233 W1123 W1223 W1233
% Rafael Henriques


Wall=zeros(3,3,3,3);

comb=[1 16 81 2 3 8 24 27 54 4 9 36 6 12 18];

for i=1:3;
    for j=1:3;
        for k=1:3;
            for l=1:3;
                Wall(i,j,k,l)=W(comb==i*j*k*l);
            end
        end
    end
end


