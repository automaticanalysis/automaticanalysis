
function F=DTImodel2(x,xdata)
%DTImodel Single diffusion tensor model for diffusion MRI analysis

v1=xdata(1,:);
v2=xdata(2,:);
v3=xdata(3,:);
v4=xdata(4,:);
v5=xdata(5,:);
v6=xdata(6,:);

%F=x(1)*exp(v1.*x(2)+v2.*x(3)+v3.*x(4)+v4.*x(5)+v5.*x(6)+v6.*x(7));
%F=x(1)*v2.*v1.^2+x(2)*sin(v1)+x(3)*v2.^3;

F=x(7)*exp(x(1)*v1/1000+x(2)*v2/1000+x(3)*v3/1000+x(4)*v4/1000+x(5)*v5/1000+x(6)*v6/1000);