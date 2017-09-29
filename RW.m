function [ V ] = RW( Vold,alpha_E,x,H,A,Psi)
%RW Rescorla-Wagner model
%   This function calculates one time-step of the Rescorla-Wagner model
%   -V is the current associative strength
%   -tau_alpha is the learning rate (inverse of alpha)
%   -z is the US trace
%   -x is the CS trace
%   -h is the time step

%Psi=max(1*10^(-6),Psi); % this takes care of very small Psi values

% if ((H*A)/Psi-max(0,x*V'))>=0
%     V=V(1)+alpha_E*((H*A)/Psi-x*V')*x(1);
% else
%     V=V(1)+alpha_I*((H*A)/Psi-x*V')*x(1);

V=Vold(1)+alpha_E*((H*A)/Psi-x*Vold')*x(1);

% if abs(V-Vold(1))>0.4
%     V=Vold(1);
% end

end

