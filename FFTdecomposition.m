%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                                 Example                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example of a quantification of slip band spacings and strains via a fast
% Fourier transform (FFT) decomposition. 
%
% Requirements:
% - Matlab R2021a
% - Image processing toolbox
% - Signal processing toolbox
% - Mtex 5.1.0 - initialise it (startup_mtex) before running this function
%
% Sections (run them sequentially):
% - PREPROCESSING
% It merges the EBSD and DIC files and builds the images of individual
% grains needed to perform the FFT decomposition. These already share the
% coordinate system, and the DIC data comes in the required format, i.e.
% images of the sample in binarised (I{1}) and greyscale (I{2}) formats
% and a variable Ith = [ethb,ethg], where ethb and ethg are the strain
% thresholds used to build I{1} and I{2}, respectively.
% 
% - FFT decomposition
% It performs the decomposition on individual grains and outputs the
% results and processed grain images.
%
% - FFT decomposition- Grain plots
% Routine to generate some plots of interest of the grains selected.
%
% - POST-ANALYSIS
% Analysis used in the original work. It optimises the orientation of the
% EBSD, calculates real spacing from the projection ones, counts slip bands
% and plots a variety of graphs.
%
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

clear
close all


%% ---------------------------- PREPROCESSING -----------------------------
%{

load('data\EBSD.mat','ebsd')                % Load EBSD data
load('data\DICimages.mat','I','Ith')        % Load DIC data: I: Images (I{1} = binarised, I{2} = greyscale)
                                            %                Ith: Image thresholds [strain]
pEBSD = 0.5;                                % EBSD pixel size [microns]
pDIC = 0.1007;                              % DIC pixel size [microns]
smth = 1;                                   % Is EBSD smoothing needed? 0 -> No, 1 -> Yes
Ath = 2;                                    % Threshold of EBSD pixels to count a grain
angth = 5;                                  % Threshold misorientation angle between two grains [deg]
wm = 5;                                     % EBSD pixels from the outer grain boundaries deleted inwards
wp = 3*pEBSD/pDIC;                          % DIC pixels from the inner grain boundaries deleted outwards
Ablth = 1000;                               % Minimum number of EBSD pixels to build reduced images
fileS = 'results_prep.mat';                 % Name of the save file

% Running preprocessing
[ebsd,grains,g,gs,gA,gsbl,gAbl,gtop,gr] = fftd_preprocessing(ebsd,I,smth,Ath,angth,wm,wp,Ablth,fileS);
svvars = {'I','Ith','ebsd','grains','smth','g','gs','gA','gsbl','gAbl','gtop','gr','pEBSD','pDIC'};
save(fileS,svvars{:},'-append');

%}


%% -------------------------- FFT decomposition ---------------------------
%{

load('results_prep.mat','gI','Ith')
spmax = 4;                              % Maximum number of active slip planes to search for
skmbl = 1;                              % Skeletonizing         Minimum branch length (to avoid dots) [pixels]
smsd = 1;                               % Smoothing             Std. dev. of 2D Gaussian smoothing kernel
dths = 0.05;                            % Finding bands in FFT  Angle resolution to search bands [deg]
thspar = ...                            %                       Parameters:
{'MinPeakHeight',4.6,...                %                       Minimum height for theta peak detection [pixel intensity]
 'MinPeakProminence',0.45};             %                       Minimum prominence (relative height) for theta peak detection [pixel intensity]
w = 3;                                  % Cropping FFT          Width of the cropped images [pixels]
sbspar = ...                            % Finding slip bands    Parameters:
{'MinPeakHeight',0.016/Ith(2),...       %                       Minimum height for peak detection [pixel intensity]
 'MinPeakProminence',0.0056/Ith(2),...  %                       Minimum prominence (relative height) for peak detection [pixel intensity]
 'MinPeakDistance',3};                  %                       Minimum distance (separation) for peak detection [pixels]
fileS = 'results_fftd.mat';             % Name of the save file

% Running FFT decomposition
[thF,thI,thSB,smeanthi,pthF,S,lambdameanp,dp,px,py,sbI,lambdap,sbImean,Ssb0,Ssbifft] =...
    fftd_decomposition(gI,spmax,skmbl,smsd,dths,thspar,w,sbspar,fileS);

%}

% Grain plots
%{

fp = 'results_prep.mat';        % Preprocessing file
fd = 'results_fftd.mat';        % FFT decomposition file
ga = 1;                         % Grains plotted (# -> grain #, [#1,#2,#3] -> grains selected, 'all' -> all grains)
dispP = {...                    % Display plots (type '...' before any figure not wanted):
    'gr'...                         % Reduced greyscale image
    'skel'...                       % Skeletonized grain
    'thsearch'...                   % Theta mean contrast to find bright bands in FFT
    'fftbands'...                   % Bright bands in FFT
    'cropped'...                    % Cropped FFTs
    'icropped'...                   % IFFTs from cropped FFTs
    'sbpeaks'...                    % Intensity spectra with slip band peaks identified
    'sbI'...                        % Slip bands identified in original image
    'sbIadj'...                     % Slip bands identified in original image (adjusted contrast)  
    'sbIl'...                       % Same, plus slip band trace at highest peak
    'sbprof'...                     % Slip band profiles comparison
    };

p = fftd_decomposition_plots(fp,fd,ga,dispP);

%}


%% ---------------------------- POST-ANALYSIS -----------------------------
%{

load('results_prep.mat','ebsd','grains','gsbl','gAbl','Ith','pDIC','pEBSD') % Load preprocessing variables 
load('results_fftd.mat','-regexp','^(?!gF$).')                              % Load all FFT decomposition variables except gF
LO = 1;                                 % Loading orientation (1 -> X, 2 -> Y, 3 -> Z)
thth = 5;                               % Threshold to flag possibly wrong angles [deg]
thflagad = [];                          % Add pairs [i,j] of slip band orientations j in grain i which are to be removed from the analysis
fileS = 'results_postanalysis.mat';     % Name of the save file
dispP = 'ols';                          % Plots (o -> orientations, l -> slip band spacings, s -> slip band strains)

% Running analysis
[thEBSD,nzEBSD,zgEBSD,mSspEBSD,nLEBSD,thdiff,meanthdiff,thflag,thdiffr,zgEBSD0] = ...     % Orientation analysis
       fftd_analysis_orientation(ebsd,grains,gsbl,LO,thSB,thth,thflagad);

for i = 1:size(thflag,1)                % Removing outliers (thflag) from each variable
    thSB(thflag(i,1),thflag(i,2)) = nan;
    thI(thflag(i,1),thflag(i,2)) = nan;
    thF(thflag(i,1),thflag(i,2)) = nan;
    pthF(thflag(i,1),thflag(i,2)) = nan;
    lambdameanp(thflag(i,1),thflag(i,2)) = nan;
    sbImean(thflag(i,1),thflag(i,2)) = nan;
    dp(thflag(i,1),thflag(i,2)) = nan;
    
    S{thflag(i,1)}(thflag(i,2),:) = nan;
    px{thflag(i,1)}{thflag(i,2)} = [];
    py{thflag(i,1)}{thflag(i,2)} = [];
    sbI{thflag(i,1)}{thflag(i,2)} = [];
    lambdap{thflag(i,1)}{thflag(i,2)} = [];
    Ssb0{thflag(i,1)}{thflag(i,2)} = [];
    Ssbifft{thflag(i,1)}{thflag(i,2)} = [];
end

e = Ith(2);                        % Converting back from greyscale to strain
[lambdamean,lambdaall,lambdameanplane,lambdastdplane,smean,sall,smeanplane,sstdplane,nsb,nlambda,ngactive] = ...     % Counts analysis
    fftd_analysis_counts(thSB,lambdameanp,lambdap,sbI,sbImean,nzEBSD,pDIC,e);

save(fileS);


%%% Plots

clrsOr = [ 0,  0.3412, 0.7176;...        % Colors for original and optimised orientations
           1,  0.8431,      0];
clrs = [   0, 0.447, 0.7410;...          % Colors for each slip plane
        0.85, 0.325, 0.0980;...
           0,   0.8,      0];

if contains(dispP,'o')      % SLIP BAND ORIENTATION PLOTS
    maxthdiff = [max(thdiffr,[],2),max(thdiff,[],2)];   % vectors [removing outliers,optimised+removing outliers]
    minthdiff = [min(thdiffr,[],2),min(thdiff,[],2)];
    MAXthdiff = maxthdiff;
    MAXthdiff(abs(minthdiff)>abs(maxthdiff)) = minthdiff(abs(minthdiff)>abs(maxthdiff));

    figure                                  % Histograms of thdiff (total)
    edgeso = -179.5:1:179.5;
    histogram(thdiffr,edgeso,'FaceColor',clrsOr(1,:))
    hold on
    histogram(thdiff,edgeso,'FaceColor',clrsOr(2,:))
    xlim([floor(min(thdiffr,[],'all')),ceil(max(thdiffr,[],'all'))])
    xlabel('$\Delta\theta [^{\circ}]$','interpreter','latex','fontsize',20)
    ylabel('Counts (slip planes)')
    legend({'Original','Optimised'},'Location','Northeast')
    set(gca,'fontsize',16,'fontname','arial')
    box on

    figure                                  % Histograms of thdiff per plane (primary, secondary, etc.)
    SP = 1:3;
    MN = [mean(abs(thdiffr(:,1:SP(end))),'omitnan');mean(abs(thdiff(:,1:SP(end))),'omitnan')]'; % Mean
    STD = [std( abs(thdiffr(:,1:SP(end))),'omitnan');std( abs(thdiff(:,1:SP(end))),'omitnan')]'; % Std. deviation
    br = bar(SP,MN);
    hold on
    X = SP'.*[1,1] + [-0.15,0.15];
    er = errorbar(X,MN,STD);
    for i = 1:2
        er(i).LineStyle = 'none';
        er(i).LineWidth = 1;
        er(i).Color = [0 0 0];  
        br(i).FaceColor = clrsOr(i,:);
        br(i).FaceAlpha = 0.6;
    end
    ylim([-1,ceil(max(MN+STD,[],'all'))])
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('Slip plane')
    ylabel('$\overline{|\Delta\theta|} [^{\circ}]$','interpreter','latex','fontsize',20)
    legend({'Original','Optimised'},'Location','Northeast')
    box on

    CS = crystalSymmetry('m-3m');           % Plots in IPF for original and optimised orientations
    sR = fundamentalSector(CS);             % The fundamental sector only
    zg0 = vector3d(zgEBSD0');               % 3d vector format
    zg0 = Miller(zg0,CS,sR);                % Miller indices format
    zg = vector3d(zgEBSD');                 % 3d vector format
    zg = Miller(zg,CS,sR);                  % Miller indices format
    msz = 8;                                % Marker size
        figure
    subplot(1,2,1)
    scatter(zg0,MAXthdiff(:,1),sR,'MarkerSize',msz)
    caxis([-8,8])
%     text(0.05,0.6,'Original','FontSize',30)
    cb = mtexColorbar;
    set(gca,'fontname','arial')
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\Delta\theta [^{\circ}]$';
    cb.Label.FontSize = 20;
    cb.Label.FontName = 'arial';
    cb.Label.Position(1) = cb.Label.Position(1) + 0.2;
        figure
    scatter(zg,MAXthdiff(:,2),sR,'MarkerSize',msz)
    caxis([-8,8])
%     text(0.05,0.6,'Optimised','FontSize',30)
    cb = mtexColorbar;
    set(gca,'fontname','arial')
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\Delta\theta [^{\circ}]$';
    cb.Label.FontSize = 20;
    cb.Label.FontName = 'arial';
    cb.Label.Position(1) = cb.Label.Position(1) + 0.2;
end


if contains(dispP,'l')      % SLIP BAND SPACINGS PLOTS
    edgesl = 0:0.5:40;

    figure                              % Histogram of mean slip band spacings per grain
    histogram(lambdamean(:,1),edgesl,'FaceColor',clrs(1,:))
    hold on
    histogram(lambdamean(:,2),edgesl,'FaceColor',clrs(2,:))
    histogram(lambdamean(:,3),edgesl,'FaceColor',clrs(3,:))
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('$\overline{\lambda}$ [$\mu$m]','interpreter','latex','fontsize',18)
    ylabel('Counts (grains)')
    legend({'Plane 1','          2','          3'},'Location','Northeast')
    xlim([edgesl(1),edgesl(end)])
    box on

    figure                              % Histogram of slip band spacings
    histogram(lambdaall{1},edgesl,'FaceColor',clrs(1,:))
    hold on
    histogram(lambdaall{2},edgesl,'FaceColor',clrs(2,:))
    histogram(lambdaall{3},edgesl,'FaceColor',clrs(3,:))
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('$\lambda$ [$\mu$m]','interpreter','latex','fontsize',18)
    ylabel('Counts (spacings)')
    legend({'Plane 1','          2','          3'},'Location','Northeast')
    box on

    figure                              % Slip band spacings vs grain area [EBSD pixels]
    scatter(gAbl*pEBSD^2/1e6,lambdamean(:,1),10,clrs(1,:),'o','filled')
    % hold on
    % scatter(gAbl*pEBSD^2/1e6,lambdamean(:,2),10,clrs(2,:),'o','filled')
    % scatter(gAbl*pEBSD^2/1e6,lambdamean(:,3),10,clrs(3,:),'o','filled')
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('Grain area [mm^2]')
    ylabel('$\overline{\lambda}$ [$\mu$m]','interpreter','latex','fontsize',18)
    box on

    figure                              % Slip band spacings vs Schmid factor [EBSD pixels]
    scatter(mSspEBSD(:,1),lambdamean(:,1),10,clrs(1,:),'o','filled')
    % hold on
    % scatter(mSspEBSD(:,2),lambdamean(:,2),10,clrs(2,:),'o','filled')
    % scatter(mSspEBSD(:,3),lambdamean(:,3),10,clrs(3,:),'o','filled')
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('Schmid factor')
    ylabel('$\overline{\lambda}$ [$\mu$m]','interpreter','latex','fontsize',18)
    box on

    figure                                  % Mean slip band spacings in IPF
    CS = crystalSymmetry('m-3m'); 
    sR = fundamentalSector(CS);             % The fundamental sector only
    nL = vector3d(nLEBSD');                 % 3d vector format
    nL = Miller(nL,CS,sR);                  % Miller indices format
    scatter(nL,lambdamean(:,1),sR,'MarkerSize',gAbl/80)
    gf = gcf;
    lambdamean1 = lambdamean(:,1);
    lambdamean1max = round(max(lambdamean1(~isinf(lambdamean1))),-1);
%     caxis([0,lambdamean1max])
    caxis([0,10])
    cb = mtexColorbar;
    set(gca,'fontname','arial')
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\overline{\lambda}$ [$\mu$m]';
    cb.Label.FontSize = 20;
    cb.Label.FontName = 'arial';
    gf.Position(3) = gf.Position(3)+20;
end

if contains(dispP,'s')      % SLIP BAND STRAINS PLOTS
%     edges = (0:1/60:1)*e;
    edgess = 0:0.004:0.25;

    figure                              % Histogram of mean slip band spacings per grain
    histogram(smean(:,1),edgess,'FaceColor',clrs(1,:))
    hold on
    histogram(smean(:,2),edgess,'FaceColor',clrs(2,:))
    histogram(smean(:,3),edgess,'FaceColor',clrs(3,:))
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('$\overline{\epsilon_{\textrm{SB}}}$ [-]','interpreter','latex','fontsize',18)
    ylabel('Counts (grains)')
    legend({'Plane 1','          2','          3'},'Location','Northeast')
    box on

    figure                              % Histogram of slip band spacings
    histogram(sall{1},edgess,'FaceColor',clrs(1,:))
    hold on
    histogram(sall{2},edgess,'FaceColor',clrs(2,:))
    histogram(sall{3},edgess,'FaceColor',clrs(3,:))
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('$\epsilon_{\textrm{SB}}$ [-]','interpreter','latex','fontsize',18)
    ylabel('Counts (slip bands)')
    legend({'Plane 1','          2','          3'},'Location','Northeast')
    box on

    figure                              % Slip band spacings vs grain area [mm2]
    scatter(gAbl*pEBSD^2/1e6,smean(:,1),10,clrs(1,:),'o','filled')
    % hold on
    % scatter(gAbl*pEBSD^2/1e6,smean(:,2),10,clrs(2,:),'o','filled')
    % scatter(gAbl*pEBSD^2/1e6,smean(:,3),10,clrs(3,:),'o','filled')
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('Grain area [mm^2]')
    ylabel('$\overline{\epsilon_{\textrm{SB}}}$ [-]','interpreter','latex','fontsize',18)
    ylim([0,edgess(end)])
    box on

    figure                              % Slip band spacings vs Schmid factor [EBSD pixels]
    scatter(mSspEBSD(:,1),sbImean(:,1)*e,10,clrs(1,:),'o','filled')
    % hold on
    % scatter(mSspEBSD(:,2),sbImean(:,2)*e,10,clrs(2,:),'o','filled')
    % scatter(mSspEBSD(:,3),sbImean(:,3)*e,10,clrs(3,:),'o','filled')
    set(gca,'fontsize',16,'fontname','arial')
    xlabel('Schmid factor')
    ylabel('$\overline{\epsilon_{\textrm{SB}}}$ [-]','interpreter','latex','fontsize',18)
    xlim([0,0.5])
    ylim([0,edgess(end)])
    box on

    figure                                  % Mean slip band spacings in IPF
    CS = crystalSymmetry('m-3m'); 
    sR = fundamentalSector(CS);             % The fundamental sector only
    nL = vector3d(nLEBSD');                 % 3d vector format
    nL = Miller(nL,CS,sR);                  % Miller indices format
    scatter(nL,smean(:,1),sR,'MarkerSize',gAbl/80)
    gf = gcf;
    caxis([0,edgess(end)])
    cb = mtexColorbar;
    set(gca,'fontname','arial')
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\overline{\epsilon_{\textrm{SB}}}$ [-]';
    cb.Label.FontSize = 20;
    cb.Label.FontName = 'arial';
    gf.Position(3) = gf.Position(3)+20;
end

if (contains(dispP,'l') && contains(dispP,'s'))      % STRAIN VS SPACING PLOTS
    figure                                  % Strain vs spacing
    errbs = {'Marker','^','MarkerSize',10,'LineWidth',1.2,'Color','k','CapSize',8};
    scatter(lambdamean(:,1),smean(:,1),7,clrs(1,:),'o','filled')
    hold on
    scatter(lambdamean(:,2),smean(:,2),7,clrs(2,:),'o','filled')
    scatter(lambdamean(:,3),smean(:,3),7,clrs(3,:),'o','filled')
    errorbar(lambdameanplane(1),smeanplane(1),sstdplane(1),sstdplane(1),lambdastdplane(1),lambdastdplane(1),'MarkerFaceColor',clrs(1,:),errbs{:})
    errorbar(lambdameanplane(2),smeanplane(2),sstdplane(2),sstdplane(2),lambdastdplane(2),lambdastdplane(2),'MarkerFaceColor',clrs(2,:),errbs{:})
    errorbar(lambdameanplane(3),smeanplane(3),sstdplane(3),sstdplane(3),lambdastdplane(3),lambdastdplane(3),'MarkerFaceColor',clrs(3,:),errbs{:})
    legend({'Plane 1','          2','          3'},'Location','Northeast')
    xlabel('$\overline{\lambda}$ [$\mu$m]','interpreter','latex','fontsize',18)
    ylabel('$\overline{\epsilon_{\textrm{SB}}}$ [-]','interpreter','latex','fontsize',18)
    ga = gca;
    ga.XLim(1) = 0;
    ylim([0,1.1* max(smean,[],'all')])
    set(gca,'fontsize',16,'fontname','arial')
    box on
end

%}

