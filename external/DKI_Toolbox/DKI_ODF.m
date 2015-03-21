function KODF=JHT(Dv,Kv,x,alfa,negate)

az=x(1);
el=x(2);



MD=mean(Dv(1:3)); % Mean diffusivity

DT=[Dv(1),Dv(4),Dv(5);...
    Dv(4),Dv(2),Dv(6);...
    Dv(5),Dv(6),Dv(3) ];

U=inv(DT)*MD;% dimensionless tensor U

[xxx, yyy, zzz]=sph2cart(az,el,1);

gnz=[xxx,yyy,zzz];

Un=U*gnz';
nUn=gnz*U*gnz';

GODF=(1/(nUn)).^((alfa+1)/2);

% matrice to compute KODF
V11= Un(1)^2/(nUn);
V22= Un(2)^2/(nUn);
V33= Un(3)^2/(nUn);
V12= Un(1)*Un(2)/(nUn);
V13= Un(1)*Un(3)/(nUn);
V23= Un(2)*Un(3)/(nUn);
%V2=[V11 V12 V13; ...
%    V12 V22 V23; ...
%    V13 V23 V33];


KODF=GODF*(1+1/24*(Kv(1)*(3*U(1,1)*U(1,1)-6*(alfa+1)*U(1,1)*V11+(alfa+1)*(alfa+3)*V11*V11)+...
    Kv(2)*(3*U(2,2)*U(2,2)-6*(alfa+1)*U(2,2)*V22+(alfa+1)*(alfa+3)*V22*V22)+...
    Kv(3)*(3*U(3,3)*U(3,3)-6*(alfa+1)*U(3,3)*V33+(alfa+1)*(alfa+3)*V33*V33)+...
    Kv(4)*(12*U(1,1)*U(1,2)-12*(alfa+1)*U(1,1)*V12-12*(alfa+1)*U(1,2)*V11+4*(alfa+1)*(alfa+3)*V11*V12)+...
    Kv(5)*(12*U(1,1)*U(1,3)-12*(alfa+1)*U(1,1)*V13-12*(alfa+1)*U(1,3)*V11+4*(alfa+1)*(alfa+3)*V11*V13)+...
    Kv(6)*(12*U(1,2)*U(2,2)-12*(alfa+1)*U(1,2)*V22-12*(alfa+1)*U(2,2)*V12+4*(alfa+1)*(alfa+3)*V12*V22)+...
    Kv(7)*(12*U(2,2)*U(2,3)-12*(alfa+1)*U(2,2)*V23-12*(alfa+1)*U(2,3)*V22+4*(alfa+1)*(alfa+3)*V22*V23)+...
    Kv(8)*(12*U(1,3)*U(3,3)-12*(alfa+1)*U(1,3)*V33-12*(alfa+1)*U(3,3)*V13+4*(alfa+1)*(alfa+3)*V13*V33)+...
    Kv(9)*(12*U(2,3)*U(3,3)-12*(alfa+1)*U(2,3)*V33-12*(alfa+1)*U(3,3)*V23+4*(alfa+1)*(alfa+3)*V23*V33)+...
    Kv(10)*(6*U(1,1)*U(2,2)+12*U(1,2)*U(1,2)-6*(alfa+1)*U(1,1)*V22-6*(alfa+1)*U(2,2)*V11-24*(alfa+1)*U(1,2)*V12+2*(alfa+1)*(alfa+3)*V11*V22+4*(alfa+1)*(alfa+3)*V12*V12)+...
    Kv(11)*(6*U(1,1)*U(3,3)+12*U(1,3)*U(1,3)-6*(alfa+1)*U(1,1)*V33-6*(alfa+1)*U(3,3)*V11-24*(alfa+1)*U(1,3)*V13+2*(alfa+1)*(alfa+3)*V11*V33+4*(alfa+1)*(alfa+3)*V13*V13)+...
    Kv(12)*(6*U(2,2)*U(3,3)+12*U(2,3)*U(2,3)-6*(alfa+1)*U(2,2)*V33-6*(alfa+1)*U(3,3)*V22-24*(alfa+1)*U(2,3)*V23+2*(alfa+1)*(alfa+3)*V22*V33+4*(alfa+1)*(alfa+3)*V23*V23)+...
    Kv(13)*(12*U(1,1)*U(2,3)+24*U(1,2)*U(1,3)-12*(alfa+1)*U(1,1)*V23-12*(alfa+1)*U(2,3)*V11-24*(alfa+1)*U(1,2)*V13-24*(alfa+1)*U(1,3)*V12+4*(alfa+1)*(alfa+3)*V11*V23+8*(alfa+1)*(alfa+3)*V12*V13)+...
    Kv(14)*(12*U(2,2)*U(1,3)+24*U(2,1)*U(2,3)-12*(alfa+1)*U(2,2)*V13-12*(alfa+1)*U(1,3)*V22-24*(alfa+1)*U(1,2)*V23-24*(alfa+1)*U(2,3)*V12+4*(alfa+1)*(alfa+3)*V22*V13+8*(alfa+1)*(alfa+3)*V12*V23)+...
    Kv(15)*(12*U(3,3)*U(1,2)+24*U(3,1)*U(3,2)-12*(alfa+1)*U(3,3)*V12-12*(alfa+1)*U(1,2)*V33-24*(alfa+1)*U(1,3)*V23-24*(alfa+1)*U(2,3)*V13+4*(alfa+1)*(alfa+3)*V33*V12+8*(alfa+1)*(alfa+3)*V13*V23)));

if nargin == 5 && negate
	KODF = -KODF;
end
%for i=1:3
%    for j=1:3
%        for k=1:3
%            for l=1:3
%                sum18=sum18+KT(i,j,k,l)*(3*U(i,j)*U(k,l)...
%                    -6*(alfa+1)*U(i,j)*V2(k,l)...
%                    +(alfa+1)*(alfa+3)*V2(i,j)*V2(k,l));
%            end
%        end
%             end
%         end
%
%         disp(toc)
%
%         KODF=GODF*(1+1/24*sum18);
%
