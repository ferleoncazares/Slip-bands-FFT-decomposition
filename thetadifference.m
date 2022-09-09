%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                Analysis                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to calculate the differences between two sets of angles (EBSD
% and DIC).
%
% Requirements:
% - Matlab R2021a
%
% Inputs:
% thEBSD        Angles of slip bands calculated from the EBSD data [deg]
% thSB          Angles of slip bands found via FFT decomposition [deg]
% thth          Threshold to flag possibly wrong angles [deg]
%
% Outputs:
% thEBSDm       thEBSD angles sorted to match those in thSB [deg]
% thdiff        Differences between slip band angles in EBSD and FFT decomposition [deg]
% meanthdiff    Mean thdiff angle difference [deg]
% thflag        [i,j] indices of grain i and trace orienation j with thdiff > thth
% ord           Order used to sort thEBSDm
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [thEBSDm,thdiff,meanthdiff,thflag,ord] = thetadifference(thEBSD,thSB,thth)

thEBSDm = zeros(size(thEBSD));
ord = zeros(size(thEBSD));
thflag = nan(size(thEBSD));

ordp = perms(1:4);
for i = 1:size(thEBSD,1)
    thEBSDp = perms(thEBSD(i,:));               % All permutations of thEBSD
    thd1 = abs(thSB(i,:).*ones(size(thEBSDp,1),1) - thEBSDp);       %... alternatives to wrap around +- 90 deg
    thd2 = abs(thSB(i,:).*ones(size(thEBSDp,1),1) - thEBSDp + 180);
    thd3 = abs(thSB(i,:).*ones(size(thEBSDp,1),1) - thEBSDp - 180);
    thdmin = min(min(thd1,thd2),thd3);
    thdp = sum(thdmin,2,'omitnan');             % th difference for all cases
    [~,thi] = min(thdp);                        % Choose the better match
    thEBSDm(i,:) = thEBSDp(thi,:);              % Build thEBSD with the best match
    ord(i,:) = ordp(thi,:);                     % Order for the matching
    thflag(i,:) = min(thdmin,[],1) > thth;      % Flagging thSB values further than thth deg from the closest thEBSD
end
thdiff = thSB - thEBSDm;                        % Angle differences
thdiff(thdiff > 90) = thdiff(thdiff > 90)-180;  %... to wrap around +-180 deg
thdiff(thdiff <-90) = thdiff(thdiff <-90)+180;

meanthdiff = mean(abs(thdiff),'all','omitnan'); % Mean angle difference

[iflag,jflag] = ind2sub(size(thflag),find(thflag));     % Rearranging thflag
[~,kflag] = sort(iflag,'ascend');
thflag = [iflag(kflag),jflag(kflag)];


end