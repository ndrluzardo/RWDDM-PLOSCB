function [ V1,alpha ] = RW_PH( V,alpha,x,alpha_E,H,A,Psi,gamma)
%PH Pearce-Hall alpha with RW error correction
%   This function calculates one time-step of the Pearce-Hall model
%   -V may be either a scalar or a row vector. If vector, place Vs in
%   order: [V1 V2 V3 ...];
%   -x is the same as V;
%   

Psi=max(0.1e-3,Psi); % this takes care of very small Psi values

error=((H*A)/Psi-x*V);

if ((H*A)/Psi-x*V)>=0
    V1=V+alpha*alpha_E*(error)*x;
else
    V1=V+alpha*alpha_E*(error)*x;
end

alpha=alpha+gamma*(abs(error)-alpha);
if alpha>1
    alpha=0.9;
end


end

