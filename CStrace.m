function [ x ] = CStrace( Input,mu,sigma,tau_x, CS, x, h )
%CStrace is a Gaussian RBF with TDDM input and a trace 
%   Detailed explanation goes here

switch CS
    case 1
        x=exp(-(Input-mu).^2/(2*sigma^2));
    case 0
        x=x+(h/tau_x)*-x;
end

