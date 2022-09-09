%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                Analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to calculate true slip band spacings and planar strains, as well
% as numbers of grains, active slip planes, slip bands and slip band
% spacings counted, and their mean values. This must be performed after
% fftd_analysis_orientation.
%
% Requirements:
% - Matlab R2021a
%
% Inputs:
% thSB          Angles of slip bands found via FFT decomposition [deg]
% lambdameanp   Mean slip band projection spacings [pixels]
% lambdap       Individual slip band projection spacings [pixels]
% sbI           Peak intensities of slip bands detected
% sbImean       Mean peak intensities of slip bands detected per grain and search line
% nzEBSD        z-component of octahedral plane unit normals in optimised EBSD
% pDIC          DIC pixel size [um]
% e             Strain to convert pixel intensities back to strains [-]
%
% Outputs:
% lambdamean        Mean slip band spacing per grain slip plane [um]
% lambdaall         All slip band spacings sorted by plane [um]
% lambdameanplane   Mean slip band spacing per plane [um]
% lambdastdplane    Std. dev. of slip band spacings per plane [um]
% smean             Mean slip band strains per grain slip plane
% sall              All slip band strains sorted by plane
% smeanplane        Mean slip band strains per plane
% sstdplane         Std. dev. of slip band strains per plane
% nsb               Numbers of slip bands detected per slip plane
% nlambda           Numbers of slip band spacings detected per slip plane
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [lambdamean,lambdaall,lambdameanplane,lambdastdplane,smean,sall,smeanplane,sstdplane,nsb,nlambda,ngactive] = ...
    fftd_analysis_counts(thSB,lambdameanp,lambdap,sbI,sbImean,nzEBSD,pDIC,e)

%% Numbers of grains and slip bands
sbIall = cell(1,4);                     % Pixel intensity of all slip bands detected (used to build nsb)
ngsb = zeros(length(sbI),4);            % Number of slip bands in each grain and slip plane
for i = 1:length(sbI)
    for j = 1:length(sbI{i})
        if ~isempty(sbI{i}{j})
            ngsb(i,j) = length(sbI{i}{j});
            sbIall{j} = [sbIall{j},sbI{i}{j}];
        end
    end
end
nsb = sum(ngsb);                        % Numbers of slip bands detected in planes [1, 2, 3, 4]
ngactive = sum(~isnan(thSB));           % Numbers of grains with slip bands detected (FFT bright bands) in planes [1, 2, 3, 4]

disp('Counts:')
disp(['   - Grains analysed: ',num2str(size(thSB,1))])
disp(['   - Grains with slip bands detected in [0,1,2,3,4] planes: [',num2str([length(lambdap)-ngactive(1),ngactive-[ngactive(2:end),0]]),']'])
disp(['   - Slip bands detected: ',num2str(sum(nsb))])
disp(['   - Slip bands detected in planes [1,2,3,4]: [',num2str(nsb),']'])


%% Slip band spacings
lambdamean = lambdameanp.*sqrt(1-nzEBSD.^2)*pDIC;                   % Mean slip band spacing per grain slip plane

lambda = cell(size(lambdap)); 
lambdaall = cell(1,4);
for i = 1:length(lambdap)
    lambda{i} = cell(size(lambdap{i}));
    for j = 1:length(lambdap{i})
        lambda{i}{j} = lambdap{i}{j}*sqrt(1-nzEBSD(i,j)^2)*pDIC;    % Slip band spacings measured in each grain and plane orientation
        lambdaall{j} = [lambdaall{j},lambda{i}{j}];                 % Same but grouped by plane order
    end
end
nlambda = cellfun(@length,lambdaall);                               % Number of slip band spacings measured per planr                         
lambdameanplane = cellfun(@mean,lambdaall);                         % Mean slip band spacing per plane...
lambdastdplane = cellfun(@std,lambdaall);                           % ... and std. dev.

disp('Slip band spacings:')
disp(['   - Slip band spacings measured: ',num2str(sum(nlambda))])
disp(['   - Slip band spacings measured in planes [1,2,3,4]: [',num2str(nlambda),']'])
disp(['   - Mean slip band spacings in planes [1,2,3,4]: [',num2str(lambdameanplane),']'])

%% Slip band strains
sall = sbIall;                                                      % All strains detected
for i = 1:length(sall)
    sall{i} = sall{i}*e;
end
smean = sbImean*e;                                                  % Mean slip band plane strain per grain per plane 
smeanplane = cellfun(@mean,sall);                                   % Mean slip band plane strain per plane...
sstdplane = cellfun(@std,sall);                                     % ... and its std. dev.
disp('Slip band planar strains:')
disp(['   - Mean slip band planar strains in planes [1,2,3,4]: [',num2str(smeanplane),']'])


end