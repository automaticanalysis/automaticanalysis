function i = checkmatlabreq(req)
i = getmatlabversion;
i(1:numel(req),2) = req;
i = conv_1000_to_dec(i);
i = i(1) >= i(2);
end

function num = conv_1000_to_dec(col)
for n = 1:size(col,2)
    thou = col(:,n);
    num(n) = 0;
    for i = 0:numel(thou)-1
        num(n) = num(n) + thou(numel(thou)-i) * 1000^i;
    end
end
end
