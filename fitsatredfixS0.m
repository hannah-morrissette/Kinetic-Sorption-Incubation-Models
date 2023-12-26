%
%  integrate saturated model
%
function yest = fitsatredfixS0(beta,tinred)
global kads kdes Qmax qm y0z modout
%
% Assign integration parameters
%
kads=beta(1);
kdes=beta(2);
Qmax = qm;
yest=zeros(27,1);
modout=zeros(482,1);

%
% set up integration for Low
% y(:,1) = Soil Carbon
% y(:,2) = Water Carbon
%
t0 = 0;
tfinal = 24;
tspan = t0:0.1:24;

%
%  y0 is the initial values of the integration S0 and W0
%
y0=y0z(:,1);
%
%  solve the ordinary differential equation calling on the satmod function
%  which specifies the ODEs
%
[t,y] = ode23(@satmod,[t0 tfinal],y0);
%
% use output to get values of W vector at input time for low
%
yest(1:12)=interp1(t,y(:,2),tinred(1:12));
modout(1:241)=interp1(t,y(:,2),tspan);

%
% set up integration for high in the same way
%
t0 = 0;
tfinal = 24;
y0=y0z(:,2);
[t,y] = ode23(@satmod,[t0 tfinal],y0);
%
% use output to get values of W vector at input time for high completing
% the set of estimates using current set of parameters
%
yest(13:27)=interp1(t,y(:,2),tinred(13:27));
modout(242:482)=interp1(t,y(:,2),tspan);

%
end