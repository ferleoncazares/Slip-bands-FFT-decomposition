%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                Analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to determine the rotation that minimises the angle differences
% between EBSD and DIC.
%
% Requirements:
% - Matlab R2021a
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Inputs:
% ebsd          EBSD dataset of surviving grains
% grains        Grains of EBSD dataset
% thSB          Angles of slip bands found via FFT decomposition [deg]
%
% Outputs:
% R             Optimal set of rotations along x-, y- and z-axes, in that order 
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function R = m_EBSD_matchingtraces(ebsd,grains,thSB)

ss = symmetrise(slipSystem.fcc(ebsd('Ni-superalloy').CS),'antipodal');  

fun = @(R) m_EBSD_thdiffR(grains,ss,thSB,R);        % Function to optimise
R0 = [0,0,0];                                       % Initial guess            

R = fminsearch(fun,R0);                             % Minimising sum(abs(thdiff))

end

function [thdiffR] = m_EBSD_thdiffR(grains,ss,thSB,R)
    
    ssgE = grains.meanOrientation * ss;                 % Slip systems for each grain
    
    Rm = rotation('axis',vector3d.Z,'angle',R(3)*degree) .* ...         % ROTATIONS in order A' = RmA = RzRxRzA
         rotation('axis',vector3d.X,'angle',R(2)*degree) .* ...
         rotation('axis',vector3d.Z,'angle',R(1)*degree);

    ssg = Rm .* ssgE;                                   % Slip systems in rotated coordinate system
    
    trEBSD = trace(ssg(:,[1,4,7,10]));                  % Traces of all slip planes
    thEBSD = atand(trEBSD(:,:).y./trEBSD(:,:).x);       % Angles of slip traces
    
    [~,~,thdiffR] = thetadifference(thEBSD,thSB,0);     % Matrix with thSB - thEBSD

%     disp(num2str(thdiffR))
end

