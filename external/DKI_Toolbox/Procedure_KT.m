function [MV,MR]=Procedure_KT(Dv,Kv,V,Nvis,Uvis,optionsT,fraqMD)
% This function finds the maxima of perpendicular kurtosis
% In white matter these directions can be used as estimates of fibre
% direction
%
% Inputs:
% Dv=[D11 D22 D33 D12 D13 D23];
% Kv=[W1111 W2222 W3333 W1112 W1113...
%     W1222 W2223 W1333 W2333 W1122...
%     W1133 W2233 W1123 W1223 W1233];
% V = Initial search grid
% Nvis = number of grid point neighbour
% Uvis = neighbour index
% optionsT = convergence options
% fraqMD = threshold to remove false fiber directions
%
% Output:
% MV = Maxima pKT directions
% MR = Maxima pKT values
%
% Developer: Rafael Neto Henriques
% First version implemented 02/07/2014
% Updates:
% O1/09/2014 - RNH - Correction not done if fibre direction estimate
%                       is only one
% 09/04/2015 - RNH - New sub-functions, Per_circle should be faster now

% step0 - Define some parametres
nodir=30;
alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
calfa=cos(alfa)';
salfa=sin(alfa)';
MD=mean(Dv(1:3));

% step 1 - Find the values for the initial search grid
R=RadialKurtosis_v(Dv,Kv,V,MD,nodir,calfa,salfa);

% step 2 - detect maxima within the initial points
np=length(R);
MeM=false(np,1);
for vt=1:np
    condi=(R(vt)>R(Uvis(vt,1:Nvis(vt))));
    if sum(condi)==Nvis(vt)
        MeM(vt)=true;
    end
end

% step 3 - select the maxima points
MV=V(MeM,:);
MR=R(MeM,:);

% step 4 - convergence
if ~isempty(MR)
    [az,el] = cart2sph(MV(:,1),MV(:,2),MV(:,3));
    lMR=length(MR);
    for vm=1:lMR
        xi(2)=el(vm);
        xi(1)=az(vm);
        myf=@(x) RadialKurtosis(Dv,Kv,x,true,MD,nodir,calfa,salfa);
        xf = fminsearch(myf,xi,optionsT);
        az(vm)=xf(1);
        el(vm)=xf(2);
        MR(vm)=RadialKurtosis(Dv,Kv,xf,false,MD,nodir,calfa,salfa);
    end
    [x,y,z] = sph2cart(az,el,1);
    MV=[x,y,z];
    
    % step 5 select only unique directions
    vaf=(sum(abs(triu(MV*MV',1))>0.99)==0);
    MV=MV(vaf,:);
    MR=MR(vaf);
    
    % step 6 correction of maxima that does not correspond to direction of
    % fibres
    np=length(MR);
    rfp=true(np,1);
    
    if length(MR)>1
        if fraqMD ~= 0
            for vm=1:length(MR)
                gn=MV(vm,:);
                DT=[Dv(1),Dv(4),Dv(5);...
                    Dv(4),Dv(2),Dv(6);...
                    Dv(5),Dv(6),Dv(3)];
                MD=mean(Dv(1:3));
                if gn*DT*gn'< fraqMD*MD;
                    rfp(vm)=false;
                end
            end
            MV=MV(rfp,:);
            MR=MR(rfp,:);
        end
    end
end