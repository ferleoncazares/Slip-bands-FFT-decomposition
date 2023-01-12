%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                Analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to calculate relevant EBSD data upon applying a rotation.
%
% Requirements:
% - Matlab R2021a
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Inputs:
% ebsd          EBSD dataset of surviving grains
% grains        Grains of EBSD dataset
% LO            Loading orientation (1 -> X, 2 -> Y, 3 -> Z)
% R             Rotations along x-, y- and z-axes in that order
% Mflag         Taylor factor is only calculated if Mflag == 1 
%
% Outputs:
% trEBSD        Traces of octahedral planes in rotated EBSD
% thEBSD        Angles of octahedral plane traces in rotated EBSD [deg]
% nzEBSD        z-component of octahedral plane unit normals in rotated EBSD
% zgEBSD        Crystallographic direction in-and-out of sample surface in rotated EBSD
% mSspEBSD      Schmid factors of octahedral planes in rotated EBSD
% MSspEBSD      Taylor factors for strain tensor diag([1,-1/3,-1/3]) in rotated EBSD
% nLEBSD        Loading direction (crystallographic orientation) in rotated EBSD
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [trEBSD,thEBSD,nzEBSD,zgEBSD,mSspEBSD,MSspEBSD,nLEBSD] = m_EBSD_gdataR(ebsd,grains,LO,R,Mflag)

ss = symmetrise(slipSystem.fcc(ebsd('Ni-superalloy').CS),'antipodal');  
Rg = grains.meanOrientation;                                        % Grain angles; 
ssgE = Rg * ss;                                                     % Slip systems for each grain

Rm = rotation('axis',vector3d.Z,'angle',R(3)*degree) .* ...         % ROTATIONS in order A' = RmA = RzRxRzA
     rotation('axis',vector3d.X,'angle',R(2)*degree) .* ...
     rotation('axis',vector3d.Z,'angle',R(1)*degree);
 
ssg = Rm .* ssgE;                                                   % Slip systems in corrected coordinate system
 
% Traces
trEBSD = trace(ssg(:,[1,4,7,10]),vector3d.Z);                       % Traces of all slip planes
thEBSD = atand(trEBSD(:,:).y./trEBSD(:,:).x);                       % Angles of slip traces
nzEBSD = abs(ssg(:,[1,4,7,10]).n.z)./ norm(ssg(:,[1,4,7,10]).n);    % z-component of slip planes

% Z-direction
zgEBSD = inv(Rg) .* inv(Rm) .* vector3d.Z;                          % z-direction of each grain
zgEBSD = sort([abs(zgEBSD.x),abs(zgEBSD.y),abs(zgEBSD.z)],2,'ascend');      %... to show in the same region of IPF
zgEBSD(:,1) = -zgEBSD(:,1);

% Schmid factors and loading orientation
if LO == 1
    Ssample = stressTensor.uniaxial(vector3d.X);    % X             % Stress in sample coordinates
    esample = strainTensor(diag([1,-1/3,-1/3]));                    % Strain tensor in sample coordinates
    nLEBSD = inv(Rg) .* inv(Rm) .* vector3d.X;                      % Loading orientation in crystal coordinates
elseif LO ==2
    Ssample = stressTensor.uniaxial(vector3d.Y);    % Y
    esample = strainTensor(diag([-1/3,1,-1/3]));
    nLEBSD = inv(Rg) .* inv(Rm) .* vector3d.Y;
elseif LO == 3
    Ssample = stressTensor.uniaxial(vector3d.Z);    % Z
    esample = strainTensor(diag([-1/3,-1/3,1]));
    nLEBSD = inv(Rg) .* inv(Rm) .* vector3d.Z;
else
    error('LO must be 1, 2 or 3!')
end
mSg = abs(ssg.SchmidFactor(Ssample,'antipodal'));                   % Schmid factors for each grain
mSspEBSD = [max(mSg(:, 1: 3),[],2),...                              % Maximum Schmid factor per slip plane per grain
            max(mSg(:, 4: 6),[],2),...
            max(mSg(:, 7: 9),[],2),...
            max(mSg(:,10:12),[],2)];
if Mflag == 1
    ssMcrystal = symmetrise(slipSystem.fcc(ebsd('Ni-superalloy').CS));  % 24 slip systems ...
    ssMsample = Rm * Rg * ssMcrystal;                               % ... in the sample coordinate system
    MSspEBSD = nan(length(grains),1);                               % Assigning vector length
    for i = 1:length(grains)
        MSspEBSD(i) = calcTaylor(esample,ssMsample(i,:));           % Taylor factors
    end
else
    MSspEBSD = [];
end

nLEBSD = sort([abs(nLEBSD.x),abs(nLEBSD.y),abs(nLEBSD.z)],2,'ascend');  %... to show in the same region of IPF
nLEBSD(:,1) = -nLEBSD(:,1);

% Ordering matrices by descending Schmid factor (EBSDu to EBSD)
[mSspEBSD,so] = sort(mSspEBSD,2,'descend');                         % Order
for k = 1:size(so,1)
    trEBSD(k,:) = trEBSD(k,so(k,:));                                    % Traces
    thEBSD(k,:) = thEBSD(k,so(k,:));                                    % Trace angles
    nzEBSD(k,:) = nzEBSD(k,so(k,:));                                    % z-component of slip planes
end

end