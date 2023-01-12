%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                Analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to optimise the EBSD orientation by minimising the difference
% between the trace angles predicted via EBSD and FFT decomposition. It
% also outputs relevant EBSD data to recover the true distances and plot
% the results.
%
% Requirements:
% - Matlab R2021a
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Inputs:
% ebsd          EBSD dataset of surviving grains
% grains        Grains of EBSD dataset
% gsbl          gs indices after removing grains with gA < Ablth
% LO            Loading orientation (1 -> X, 2 -> Y, 3 -> Z)
% thSB          Angles of slip bands found via FFT decomposition [deg]
% thth          Threshold to flag possibly wrong angles [deg]
% thflagad      Add pairs [i,j] of slip band orientations j in grain i which are to be removed from the analysis
%
% Outputs:
% thEBSD        Angles of octahedral plane traces in optimised EBSD [deg]
% nzEBSD        z-component of octahedral plane unit normals in optimised EBSD
% zgEBSD        Crystallographic direction in-and-out of sample surface in optimised EBSD
% mSspEBSD      Schmid factors of octahedral planes in optimised EBSD
% MSspEBSD      Taylor factors for strain tensor diag([1,-1/3,-1/3]) in rotated EBSD
% nLEBSD        Loading direction (crystallographic orientation) in optimised EBSD
% thdiff        Differences between slip band angles in optimised EBSD and FFT decomposition
% meanthdiff    Mean thdiff angle difference
% thflag        [i,j] indices of grain i and trace orienation j with thdiff > thth
% thdiffr       Differences between slip band angles in ORIGINAL EBSD and FFT decomposition
% zgEBSD0       Crystallographic direction in-and-out of sample surface in ORIGINAL EBSD
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [thEBSD,nzEBSD,zgEBSD,mSspEBSD,MSspEBSD,nLEBSD,thdiff,meanthdiff,thflag,thdiffr,zgEBSD0] = ...
         fftd_analysis_orientation(ebsd,grains,gsbl,LO,thSB,thth,thflagad)

disp('Orientation analysis...')
ticor = tic;                                                   % Function timer

grainsbl = grains(gsbl);                % Data for surviving grains
ebsdbl = ebsd(grainsbl);

% Before optimising
[~,thEBSD0,~,zgEBSD0,~,~,~] = m_EBSD_gdataR(ebsdbl,grainsbl,LO,[0,0,0],0);  % Traces, z-vector and Schmid factors
[~,~,meanthdiff0,~] = thetadifference(thEBSD0,thSB,thth);                   % Differences in DIC and not optimised EBSD angles

% Optimised
Ro = m_EBSD_matchingtraces(ebsdbl,grainsbl,thSB);                           % Optimal rotation operations  
[~,thEBSDo,~,~,~,~,~] = m_EBSD_gdataR(ebsdbl,grainsbl,LO,Ro,0);             % Repeat operations
[~,~,meanthdiffo,thflag,~] = thetadifference(thEBSDo,thSB,thth);
thflag = [thflag;thflagad];

% Outliers removed
thSBo = thSB;
for i = 1:size(thflag,1)                % Removing outliers and rearranging thSB
    thSBo(thflag(i,1),thflag(i,2)) = nan;
end
[~,thdiffr,meanthdiffr,~,~] = thetadifference(thEBSD0,thSBo,thth);                  % Differences in DIC and not optimised EBSD angles

% Optimised + outliers removed
R = m_EBSD_matchingtraces(ebsdbl,grainsbl,thSBo);                                   % Optimal rotation operations
[trEBSD,thEBSD,nzEBSD,zgEBSD,mSspEBSD,MSspEBSD,nLEBSD] = m_EBSD_gdataR(ebsdbl,grainsbl,LO,R,1);  % Repeat operations and calculate Taylor factors
[thEBSD,thdiff,meanthdiff,~,ord] = thetadifference(thEBSD,thSBo,thth);
for i = 1:size(thdiff,1)                                                            % Rearranging to match the order in thSB
    trEBSD(i,:) = trEBSD(i,ord(i,:));
    nzEBSD(i,:) = nzEBSD(i,ord(i,:));
    mSspEBSD(i,:) = mSspEBSD(i,ord(i,:));
end

tocor = toc(ticor);
disp(['... done! Time: ',num2str(tocor),' s.'])

disp(['   - Optimal set of rotations: [',num2str(R),'] deg.'])
disp('   - Mean theta-difference values:')
disp(['     - Before optimising:      ',num2str(meanthdiff0),'  deg'])
disp(['     - Optimised:              ',num2str(meanthdiffo),' deg'])
disp(['     - ',num2str(size(thflag,1)),' outliers removed:     ',num2str(meanthdiffr),' deg'])
disp(['     - Optimised - ',num2str(size(thflag,1)),' outliers: ',num2str(meanthdiff),' deg'])

end