function M = reshapeWholeBrain(ss,D)

if numel(size(D))==2
    M = reshape(D', ss(1), ss(2), ss(3), ss(4)); 
end

if numel(size(D))==4
     M = double(reshape(D, prod(ss(1:3)),ss(4))');
end
            
            
           