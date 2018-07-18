function Pdir=Perp_circle_x(a,b,c,nodir,calfa,salfa)
% Rafael Henriques - Compute perdendicular directions relative to a vector
% of coordinates a,b, and c

if nargin<5
alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
calfa=cos(alfa)';
salfa=sin(alfa)';
end 

sq=sqrt(b^2+c^2);
Pdir=[-sq*salfa, (-c*calfa+b*a*salfa)/sq, (b*calfa+c*a*salfa)/sq];
% R1=[0, -c/sq, b/sq];
% R2=[-sq,b*a/sq,c*a/sq];
% Dir_samples=zeros(nodir,3);
% for s=1:nodir
%     Dir_samples(s,:)=R1*(calfa(s))+R2*(salfa(s));
% end
%hold all
%plot3(Dir_samples(:,1),Dir_samples(:,2),Dir_samples(:,3))
 


 


