%% 5B3D-Finger-Model-Results: cartForces
%% Author: Alexander B. Ambrose
%% Date: 8-6-2025

%% Description:
%  Assembles the cartesian grasp forces given a reference total enveloping
%  force in N and cartesian force vectors from each pad
%  Inputs:
%  F -> the vector of total enveloping forces
%  Nx -> vector of the x-direction contact forces
%  Ny -> vector of the x-direction contact forces
%  ref -> the reference forces you are testing
%
%  Outputs:
%  Fx -> scalar of the x-direction contact force for each reference force
%  Fy -> vector of the y-direction contact force for each reference force



function [fx, fy] = catForces(F, Nx, Ny, ref)
for ii = 1:length(ref)
    % find finger states that produce force close to the reference
    f = find(abs(F-ref(ii)) < 1);
    idx(:, ii) = f(1:4); 
    % linearly interpolate
    E(1, ii)= ref(ii)-F(idx(1, ii));
    E(2, ii)= F(idx(end, ii))-ref(ii);
    d = E(1, ii)/(E(1, ii)+E(2, ii));
    % approximate resulting cartesian forces
    fx(ii) = Nx(idx(1, ii)) + d*(Nx(idx(end, ii))-Nx(idx(1, ii)));
    fy(ii) = Ny(idx(1, ii)) + d*(Ny(idx(end, ii))-Ny(idx(1, ii)));
end