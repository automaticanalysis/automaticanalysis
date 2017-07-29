function Pdir=Perp_circle_y(a,b,c,nodir,calfa,salfa)
% Rafael Henriques - Compute perdendicular directions relative to a vector
% of coordinates a,b, and c

if nargin<5
alfa=pi/nodir:pi/nodir:pi;% half of points because is symetric
calfa=cos(alfa)';
salfa=sin(alfa)';
end 

sq=sqrt(a^2+c^2);
Pdir=[-(c*calfa +a*b*salfa)/sq, salfa*sq, (a*calfa-c*b*salfa)/sq];
% R1=[-c/sq, 0, a/sq];
% R2=[-a*b/sq,sq,-c*b/sq];
% Dir_samples=zeros(nodir,3);
% for s=1:nodir
%     Dir_samples(s,:)=R1*(calfa(s))+R2*(salfa(s));
% end
%hold all
%plot3(Dir_samples(:,1),Dir_samples(:,2),Dir_samples(:,3))
 


