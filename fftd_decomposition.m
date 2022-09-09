%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                              Decomposition                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to quantify slip band spacings and strains via a fast Fourier
% transform (FFT) decomposition. Identification of active slip planes is
% performed on binarised images and slip band detection on greyscale
% images. The function saves the reduced images one by one in a specified
% file to avoid potential memory issues from loading them all at the same
% time.
%
% Requirements:
% - Matlab R2021a
% - Image processing toolbox
% - Signal processing toolbox
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Inputs:
% gI        Images saved onto fileS from the preprocessing step
%  _.r      Reduced images of the grain
%  _.out    Mask of pixels within gI{i}.r outside of grain i
%  _.g      Image of gI{i}.r after applying the mask gI{i}.out
% spmax     Maximum number of active slip planes to search for
% skmbl     Skeletonizing - Minimum branch length (to avoid dots) [pixels]
% smsd      Smoothing - Std. dev. of 2D Gaussian smoothing kernel
% dths      Finding bands in FFT - Angle resolution to search bands [deg]
% thspar    Finding bands in FFT - findpeaks parameters (MinPeakHeight, MinPeakProminence, etc.)
% w         Cropping FFT - Width of the cropped images [pixels]
% sbspar    Finding slip bands - findpeaks parameters (MinPeakHeight, MinPeakProminence, etc.)
% fileS     Name of the save file
%
% Outputs:
% thF       Angles of the bright bands in the FFT [deg]
% thI       Angles of search lines in IFFT [deg]
% thSB      Angles of slip bands [deg]
% smeanthi  Spectrum of mean contrast in FFT as a function of the search angle
% pthF      Peak intensities of bright bands in FFT
% S         Intensity spectrum along search lines in IFFT
% lambdameanp   Mean slip band projection spacings [pixels]
% dp        Number of pixels within grain i along sampling line j
% pX        x-coordinates of slip bands along search line [pixels]
% pY        y-coordinates of slip bands along search line [pixels]
% sbI       Peak intensities of slip bands detected
% lambdap   Individual slip band projection spacings [pixels]
% sbImean   Mean peak intensities of slip bands detected per grain and search line
% Ssb0      Slip band intensity profiles interpolated from original image
% Ssbifft   Slip band intensity profiles interpolated from IFFT
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [thF,thI,thSB,smeanthi,pthF,S,lambdameanp,dp,pX,pY,sbI,lambdap,sbImean,Ssb0,Ssbifft] = ...
         fftd_decomposition(gI,spmax,skmbl,smsd,dths,thspar,w,sbspar,fileS)

%% Setting timer and variables
disp('FFT decomposition...')
ticfft = tic;                                                   % Function timer

if ~isempty(fileS)
    gF = cell(size(gI));                                        % Cell with images related to the FFT decomposition
    save(fileS,'gF','-v7.3')                                    % Output file
    f = matfile(fileS,'Writable',true);                         % To save partial variables (reducing memory used)
end
[thF,pthF,thI,thSB,lambdameanp,dp,sbImean] = deal(nan(length(gI),spmax));    % Non-image outputs
[smeanthi,S,pX,pY,sbI,lambdap,Ssb0,Ssbifft] = deal(cell(size(gI)));
warning('off','signal:findpeaks:largeMinPeakHeight')

%% Grain loop
wb = waitbar(0,'FFT decomposition');
for i = 1:length(gI)
    if sum(sum(gI{i}.g(:,:,1))) == 0                            % In case the binary image is empty
        warning(['Empty binary image for grain ',num2str(i),'! Skipping this grain.'])
        continue
    end
    sz = size(gI{i}.g(:,:,1));                                  % Size of i images
    [y0,x0] = deal(ceil((sz(1)+0.1)/2),ceil((sz(2)+0.1)/2));    % Pixel at the centre of the FFT

    
    %% FFT
    gFi.If = bwskel(logical(gI{i}.g(:,:,1)),'MinBranchLength',skmbl);   % Skeletonized - m.b.l. to avoid dots
    Ifftds = abs(log2(fftshift(fft2(gFi.If))));                 % FFT (from binarised) where band search is performed (real)
    Ifftds(Ifftds>Ifftds(y0,x0)) = Ifftds(y0,x0);               % Removing high contrast pixels for better visualisation
    gFi.Ifftd = Ifftds;                                         % FFT - for display only
    gFi.Ifftd = gFi.Ifftd/gFi.Ifftd(y0,x0);                     % Normalising Ifftd for better visualisation
    if smsd > 0
        Ifftds = imgaussfilt(Ifftds,smsd,'FilterDomain','spatial'); % 2D Gaussian smoothing kernel with smsd std. dev.
    end 
    gFi.Ifft = fftshift(fft2(gI{i}.g(:,:,2)));                  % FFT (complex, from greyscale)
    
    
    %% Finding bright bands in FFT
    ths = (0:dths:(180-dths));                                  % Vector with search angles
    rs = -(floor(min(sz))/2-1):(floor(min(sz))/2-1);            % Vector with search line radii
    [Rs,THs] = meshgrid(rs,ths);                                % Map of rs and ths (polar coordinates)
    Xs = Rs.*cosd(THs) + x0;                                    % x and y coordinates of such map
    Ys = Rs.*sind(THs) + y0;
    [X,Y] = meshgrid(1:sz(2),1:sz(1));                          % x and y coordinates of image pixels
    sth = interp2(X,Y,Ifftds,Xs,Ys);                            % Interpolating at points in polar coordinates
    smeanthi{i} = mean(sth,2);                                  % Mean radial contrast in FFT for varied angles
    [pthFs,thFs] = findpeaks(...                                % Find peaks (pthF - peak intensity, thi - band angles in FFT)...
        [smeanthi{i};smeanthi{i}],[ths,(180+ths)],...           % ... arranged over twice the search space so it can find peaks near ths = 0 = 180 deg
        thspar{:},...                                           % ... using the search parameters defined by the user
        'SortStr','descend');                                   % ... and sorting the FFT bands by brightness
    if isempty(thFs)
        warning(['No lines detected in grain ',num2str(i),'! Skipping this grain.'])
        continue
    end
    thFsu = mod(thFs(thFs>=90 & thFs<=(270-dths)),180);         % To include every peak only once
    pthFsu = pthFs(thFs>=90 & thFs<=(270-dths))';
    if length(thFsu) > spmax
        warning(['FFT of grain ',num2str(i),' shows ',num2str(length(thFsu)),' bright bands! Considering only ',num2str(spmax),'.'])
        thFsu = thFsu(1:spmax);
        pthFsu = pthFsu(1:spmax);
    end
    nthi = length(thFsu);                                       % Number of bright bands found
    thF(i,:) = [thFsu,nan(1,spmax-nthi)];                       % Angles of the bright bands in the FFT [deg]  
    pthF(i,:) = [pthFsu,nan(1,spmax-nthi)];
    
    mi = tand(thF(i,:));                                        % Slopes of lines in FFT
    m = mi*sz(2)/sz(1);                                         % Slopes of lines in original image
    thI(i,:) = atand(m);                                        % Angles of search lines in IFFT [deg]
    thSB(i,:) = atand(tand(thI(i,:)+90));                       % Angles of slip bands (tan & atan to take it between -90 and 90 deg)

    
    %% Cropping FFT and performing IFFT
    [gFi.Ifftdc,gFi.Iifftc] = deal(zeros(sz(1),sz(2),nthi));
    for j = 1:nthi
        dl = abs(mi(j)*X-Y+(y0-mi(j)*x0))./sqrt(mi(j)^2+(-1)^2);% Distance [pixels] to line from each pixel
        F = gFi.Ifftd;                                          % FFT for display (real)
        F(dl>w/2) = 0;                                          % Cropping
        gFi.Ifftdc(:,:,j) = F;                                  % Assigning cropped
        F = gFi.Ifft;                                           % FFT (complex)
        F(dl>w/2) = 0;                                          % Cropping
        gFi.Iifftc(:,:,j) = abs(ifft2(F)).*~gI{i}.out;          % IFFTs of the cropped FFTs (masked to show contrast only within grain i)
    end
    
    
    %% Finding slip bands
    rp = -(ceil(norm(sz))/2):(ceil(norm(sz))/2);                % Line across images
    xp = rp.*cosd(thI(i,:)') + x0;                              % x and y coordinates of lines in IFFT
    yp = rp.*sind(thI(i,:)') + y0;
    S{i} = nan(size(yp));
    [pX{i},pY{i},sbI{i},lambdap{i}] = deal(cell(1,nthi));
    for j = 1:nthi
        S{i}(j,:) = interp2(X,Y,gFi.Iifftc(:,:,j),xp(j,:),yp(j,:));    % Intensity spectrum - interpolating at points in the sampling line
        [sbI{i}{j},prp] = findpeaks(S{i}(j,:),...               % Find peaks (sbI - peak intensity, prp - location along rp vector)
            sbspar{:});                                         % ... using the search parameters defined by the user
        pX{i}{j} = xp(j,prp);                                   % Peak x-coordinates
        pY{i}{j} = yp(j,prp);                                   % Peak y-coordinates
        dp(i,j) = sum(~(isnan(S{i}(j,:)) | S{i}(j,:)==0));      % Number of pixels within grain i along sampling line j
        lambdameanp(i,j) = dp(i,j)/length(pX{i}{j});            % Mean slip band projection spacing [pixels]
        
        prnan = rp(find(diff(S{i}(j,:)==0 | isnan(S{i}(j,:)))==1) + 1); % Location of segments outside of the grain of interest (1 point per segment)
        [pr,knan] = sort([prnan,rp(prp)],'ascend');             % Combined locations of peaks and zeros
        pr(knan<=length(prnan)) = nan;                          % Turn zeros into nans
        lambdap{i}{j} = pr(2:end)-pr(1:end-1);                  % Slip band projection spacings [pixels]
        lambdap{i}{j}(isnan(lambdap{i}{j})) = [];               % Removing nans
    end
    
    
    %% Slip band deformation profiles
    [Ssb0{i},Ssbifft{i}] = deal(cell(size(sbI{i})));
    for j = 1:length(sbI{i})
        if ~isempty(sbI{i}{j})
            sbImean(i,j) = mean(sbI{i}{j});
            xk = rp.*cosd(thSB(i,j)) + pX{i}{j}';               % Sampling lines referenced to each slip band
            yk = rp.*sind(thSB(i,j)) + pY{i}{j}';
            Ssb0{i}{j}    = interp2(X,Y,gI{i}.g(:,:,2),xk,yk);  % Intensity interpolation from original image 
            Ssbifft{i}{j} = interp2(X,Y,gFi.Iifftc(:,:,j),xk,yk);   % Intensity interpolation from IFFT
        end
    end     
    
    
    %% Saving images into gF vector
    if ~isempty(fileS)
        f.gF(i,1) = {gFi};
    end
    
    waitbar(i/length(gI),wb);
end
warning('on','signal:findpeaks:largeMinPeakHeight')
close(wb)
tocfft = toc(ticfft);
disp(['... done! Grains: ',num2str(length(gI)),'; Time: ',num2str(tocfft),' s.'])


%% Saving non-image variables
svvars = {'spmax','skmbl','smsd','dths','thspar','w','sbspar',...                                                            % Inputs
          'thI','thF','thSB','smeanthi','pthF','S','lambdameanp','dp','pX','pY','sbI','lambdap','sbImean','Ssb0','Ssbifft'};   % Outputs
if ~isempty(fileS)
    save(fileS,svvars{:},'-append');
end

end
