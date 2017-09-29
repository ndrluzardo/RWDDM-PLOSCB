function [ P ] = DDM( P, A, h, m, N )
%DDM Drift-Diffusion Model
%   Calculates one time-step of the DDM
P=(P+A*h+m*sqrt(A*h)*N)*((P+A*h+m*sqrt(A*h)*N)>=0);
end

