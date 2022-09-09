%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                              Preprocessing                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to build images of individual grains needed for the FFT
% decomposition analysis. This function takes as an input an ebsd dataset
% already aligned with the pixels of the micrograph from which slip bands
% will be detected. The function saves the reduced images one by one in a
% specified file to avoid potential memory issues from loading them all at
% the same time.
%
% Requirements:
% - Matlab R2021a
% - Image processing toolbox
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Inputs:
% ebsd      EBSD dataset already aligned with I
% I         Images from which slip bands are detected ({1} -> binarised and {2} -> greyscale)
% smth      Is EBSD smoothing needed? 0 -> No, 1 -> Yes
% Ath       Threshold to count a grain [EBSD pixels]
% angth     Threshold of misorientation angle between two grains [deg]
% wm        Pixels from the outer grain boundaries deleted inwards [EBSD pixels]
% wp        Pixels from the inner grain boundaries deleted outwards [EBSD pixels]
% Ablth     Minimum number of pixels to build reduced images [EBSD pixels]
% fileS     Name of the save file
%
% Outputs:
% ebsd      EBSD dataset denoised, smoothed and fitted to the size of I
% grains    Grains in ebsd
% g         EBSD datasets of individual grains in a global coordinate system
% gs        Indices of grains sorted by size (e.g. gs(1) = largest grain)
% gA        Areas of grains sorted by size [EBSD pixels]
% gsbl      gs indices after removing grains with gA < Ablth
% gAbl      gA areas after removing grains with gA < Ablth
% gtop      Topology vector. It indicates which grains are inside of each grain using gs indices
% gr        EBSD datasets of individual grains in a local coordinate system
% gI        Images saved onto fileS
%  _.r      Reduced images of the grain
%  _.out    Mask of pixels within gI{i}.r outside of grain i
%  _.g      Image of gI{i}.r after applying the mask gI{i}.out
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function [ebsd,grains,g,gs,gA,gsbl,gAbl,gtop,gr] = fftd_preprocessing(ebsd,I,smth,Ath,angth,wm,wp,Ablth,fileS)

%% Start timer and detect input errors
tprep = tic;                                % Function timer

for i = 1:(length(I)-1)                     % Making sure all images have the same size
    if ~isequal(size(I{i}),size(I{end}))
        error('All images "I" must have the same size!')
    end
end
if ~(smth == 0 || smth == 1)
    error('smth must have a value of 0 or 1')
end

%% Matching EBSD-DIC and smoothing
if smth == 1
    disp('Matching and smoothing EBSD...')
else
    disp('Matching EBSD...')
end
ebsd = ebsd('indexed');                     % Remove unindexed
ebsd(ebsd.x < 1) = [];                      % EBSD cropped to match DIC
ebsd(ebsd.x > size(I{1},2)) = [];
ebsd(ebsd.y < 1) = [];
ebsd(ebsd.y > size(I{1},1)) = [];
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',angth*degree);  % Calculate grains

if smth == 1                                % EBSD smoothing = YES
    ebsd = smooth(ebsd,splineFilter,'fill',grains);             % Smoothing filter 
    toc1 = toc(tprep);
else                                        % EBSD smoothing = NO
    toc1 = toc(tprep);
end
disp(['... done! Grains: ',num2str(length(grains)),'; Time: ',num2str(toc1),' s.'])


%% g-vector contaning separate grains
disp('Creating g-vector...')
[cts,gId] = groupcounts(ebsd.grainId);      % Counts cts of EBSD pixels with grainId gId
cts(isnan(gId)) = [];                       % Remove NaNs
gId(isnan(gId)) = [];
gId(cts<Ath) = [];                          % Remove grains with less pixels than ctsth
cts(cts<Ath) = [];
[gA,gi] = sort(cts,'descend');              % Sorted number of EBSD pixels per grain
gs = gId(gi);                               % Index of grains ordered by size (e.g. gs(1) = largest grain)

g = cell(size(gs));                         % Vector with the EBSD data of surviving grains
wb = waitbar(0,'2/5 Creating g-vector');
for i = 1:length(gs)
    g{i} = ebsd(ebsd.grainId == gs(i));     % EBSD data of grain i
    if mod(i,20)==0, waitbar(i/length(gs),wb); end
end
close(wb)
toc2 = toc(tprep);
disp(['... done! Grains: ',num2str(length(g)),'; Time: ',num2str(toc2-toc1),' s.'])


%% Grain topology
disp('Calculating grain topology...')
gtop = cell(size(g));                       % Topology vector: it indicates which grains are inside of each grain
[boxx,boxy] = deal(zeros(length(g),2));   
for i = 1:length(g)                         % Box (x and y) that contains each grain
    boxx(i,:) = [min(g{i}.x),max(g{i}.x)];
    boxy(i,:) = [min(g{i}.y),max(g{i}.y)];
end
wb = waitbar(0,'3/5 Calculating grain topology');
for i = 1:length(g)                         % Grain topology
    gsearch = find((boxx(i+1:end,1) > boxx(i,1)) & ...          % Which grains are inside of box i?
                   (boxx(i+1:end,2) < boxx(i,2)) & ...
                   (boxy(i+1:end,1) > boxy(i,1)) & ...
                   (boxy(i+1:end,2) < boxy(i,2)))+i;
    if ~isempty(gsearch)
        b = boundary(g{i}.x,g{i}.y,1);                          % Boundary of grain i
        gj = zeros(length(gsearch),2);                          % First pixel for each grain in box i
        for j = 1:length(gsearch)      
            gj(j,:) = [g{gsearch(j)}.x(ceil(end/2)),...         % Using pixel ~(end/2) of j grains
                       g{gsearch(j)}.y(ceil(end/2))];  
        end
        jin = inpolygon(gj(:,1),gj(:,2),g{i}.x(b),g{i}.y(b));   % Are grains gsearch within actual boundary of grain i?
        gtop{i} = gsearch(jin);                                 % Grains in gs within grain i
    end
    if mod(i,20)==0, waitbar(i/length(g),wb); end
end
close(wb)
toc3 = toc(tprep);
disp(['... done! Grains: ',num2str(length(g)),'; Time: ',num2str(toc3-toc2),' s.'])


%% Removing grain boundary pixels
disp('Removing boundary pixels...')
gss = gs;
g(gA < Ablth) = [];                         % Removing grains that do not pass the threshold
gss(gA < Ablth) = [];
gtop(gA < Ablth) = [];
gl = length(g);                             % Number of grains on which boundary removal is performed

gsbl = nan(size(g));                        % Surviving grains
gAbl  = nan(size(g));
wb = waitbar(0,'4/5 Removing boundary pixels');
for i = 1:length(gss)                       % Loop for each grain
    for j = 1:wm                            % Loop to eliminate pd pixels from the border
        k = boundary(g{i}.x,g{i}.y,1);      % Indices of points that correspond to the border
        g{i}(k) = [];                       % Discard the border
    end
    if length(g{i}) >= Ablth                % Grains above the threshold remaining pixels                     
        gsbl(i) = gss(i);                   % Index of grains ordered by size
        gAbl(i) = length(g{i});             % Surviving number of pixels
    end
    if mod(i,5)==0, waitbar(i/length(gss),wb); end
end

% Deleting grains (from all arrays) that did not survive the border reduction
g(isnan(gsbl(:)),:) = [];
gtop(isnan(gsbl(:)),:) = [];
gsbl(isnan(gsbl(:)),:) = [];
gAbl(isnan(gAbl(:)),:) = [];       
% ebsd_bl = cat(1,g{:});                    % Borderles EBSD
toc4 = toc(tprep);
close(wb)
disp(['... done! Grains: ',num2str(gl),'; Time: ',num2str(toc4-toc3),' s.'])


%% Creating reduced images
disp('Creating reduced images...')
gI = cell(size(g,1),1);                                         % Cell with grain images
save(fileS,'gI','-v7.3')                                        % Output file
f = matfile(fileS,'Writable',true);                             % To save partial variables (reducing memory used)

gr = g;                                                         % g-vector with local coordinate system
wb = waitbar(0,'5/5 Creating reduced images');
for i = 1:length(gr)
    b = boundary(gr{i}.x,gr{i}.y,1);                            % Boundary of grain pixels
    boxx = [floor(min(gr{i}.x)),ceil(max(gr{i}.x))];            % Box that contains the grain [x,y]
    boxy = [floor(min(gr{i}.y)),ceil(max(gr{i}.y))];
    [gIi.r,gIi.g] = deal(zeros(boxy(2)-boxy(1)+1,boxx(2)-boxx(1)+1,length(I)));  
    for j = 1:length(I)
        gIi.r(:,:,j) = I{j}((boxy(1):boxy(2)),(boxx(1):boxx(2)));   % Reduced images
    end
    Irx = zeros(1,size(gIi.r(:,:,1),1))+(1:size(gIi.r(:,:,1),2))';  % Pixel coordinates of reduced image [xr,yr]
    Iry = zeros(size(gIi.r(:,:,1),2),1)+(1:size(gIi.r(:,:,1),1)); 
    gr{i}.x = gr{i}.x-boxx(1);                                  % EBSD pixel coordinates of grain in reduced image [xr,yr]
    gr{i}.y = gr{i}.y-boxy(1);                         
    iin = inpolygon(Irx,Iry,gr{i}.x(b),gr{i}.y(b));             % Pixels within GB

    iout = zeros(size(Irx));      
    for j = 1:length(gtop{i})                                   % Pixels from grains within grain i
        goutx = ebsd(ebsd.grainId == gs(gtop{i}(j))).x-boxx(1);
        gouty = ebsd(ebsd.grainId == gs(gtop{i}(j))).y-boxy(1);
        bout = boundary(goutx,gouty,1);                         % Boundary of grain  within grain i
        if ~isempty(bout)
            if wp > 0                           % Additional layer removed
                po = polybuffer([goutx(bout),gouty(bout)],'lines',wp,'JointType','round');
                bfseg = [0;find(isnan(po.Vertices(:,1)));size(po.Vertices,1)+1];    % Segments created by polybuffer (e.g. (bfseg(1)+1):(bfseg(2)-1))
                [~,bfout] = min(po.Vertices(:,1));                                  % Sequence to find the one that corresponds to the outer boundary
                [~,bfo] = sort([bfout;bfseg]);                                      % ... i.e. bfo(1)
                polj = [po.Vertices((bfseg(bfo(1)-1)+1):(bfseg(bfo(1))-1),:);...
                        po.Vertices((bfseg(bfo(1)-1)+1),:)];
            else                                % No additional layer removed
                polj = [goutx(bout),gouty(bout)];
            end
            ioutj = inpolygon(Irx,Iry,polj(:,1),polj(:,2));     % DIC pixels within grain j
            iout = iout | ioutj;                                % Update iout
        end
    end
    
    gIi.out = ~iin' | iout';                                    % Pixels that do not belong to grain i
    for j = 1:length(I)
        gIi.g(:,:,j) = gIi.r(:,:,j);                            % Reduced images of grain i ONLY
        gIi.g(:,:,j) = gIi.g(:,:,j).*~gIi.out;                  % Remove unwanted pixels
    end
    
    f.gI(i,1) = {gIi};                                          % Saving images for grain i
    waitbar(i/length(gr),wb);
end
close(wb)
toc5 = toc(tprep);
disp(['... done! Grains: ',num2str(length(g)),'; Time: ',num2str(toc5-toc4),' s.'])

%% End
disp(['Preprocessing done! Grains: ',num2str(length(g)),'; Time: ',num2str(toc5),' s.'])

end