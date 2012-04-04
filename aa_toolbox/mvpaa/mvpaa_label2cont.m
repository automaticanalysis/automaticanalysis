function cont = mvpaa_label2cont(label, mode)
% Label = the labels of the different conditions
% Mode = are the labels discrete (0) or continuous (1) or binary (2)
% NOTE:
% binary is a special type, which tests if within category correlations for
% 1 are greater than correlations for 0... within categories for category 0
% are NaNed..

X = repmat(label', [1 length(label)]);
Y = repmat(label, [length(label) 1]);

cont = X - Y;
cont = abs(cont);

if mode == 0 || mode == 2
    % If discrete, only the same label gets high similarity, else low
    cont(~isnan(cont)) = ~cont(~isnan(cont));
    cont = 2*cont - 1;
elseif mode == 1
    % If continuous, the closer the value to the label, the more similar
    cont = -cont;
    cont = cont - min(cont(:))/2;
end
if mode == 2
    keepM = zeros(size(cont));
    keepM(cont<0) = 1;
    keepM(logical(~label), logical(~label)) = 1;
    cont(logical(~keepM)) = NaN;
end