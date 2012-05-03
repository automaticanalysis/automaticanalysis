function persisText(txt,ii)
% CUMDISP persistent disp
% cumdisp; initializes persistent display
% cumdisp(text); displays persistent text
%

if nargin==2
    if ii == 1
        fprintf('\n');
    end
end

persistent oldtxt;
if nargin<1
    oldtxt=''; 
    fprintf(1,'\n'); 
else
    fprintf(1,[repmat('\b',[1,length(oldtxt)]),txt]);
    oldtxt=sprintf(txt);
end