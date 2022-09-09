%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        FFT decomposition method                         % 
%                           Decomposition plots                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to plot images from the results of fftd_decomposition.
%
% Requirements:
% - Matlab R2021a
%
% Inputs:
% fp        Preprocessing file
% fd        FFT decomposition file
% ga        Grains plotted (# -> grain #, [#1,#2,#3] -> grains selected, 'all' -> all grains)
% dispP     Plots to be displayed.
%           Input as: dispP = {'skel','thsearch','fftbands'};  
%
% Plots:
% gr        Reduced greyscale image
% skel      Skeletonized grain
% thsearch  Theta mean contrast to find bright bands in FFT
% fftbands  Bright bands in FFT
% cropped   Cropped FFTs
% icropped  IFFTs from cropped FFTs
% sbpeaks   Intensity spectra with slip band peaks identified
% sbI       Slip bands identified in original image
% sbIadj    Slip bands identified in original image (contrast adjusted)
% sbIl      Same, plus slip band trace at highest peak
% sbprof    Slip band profiles comparison
%
% Output:
% p         Cell with all the plots generated
%
% Coded by F.D. León-Cázares
% https://orcid.org/0000-0002-3828-6695
% https://www.researchgate.net/profile/Fernando-Daniel-Leon-Cazares
%

function p = fftd_decomposition_plots(fp,fd,ga,dispP)

load(fd,'-regexp','^(?!gF$).')      % Loading results from fft decomposition
fp = matfile(fp);                   % To load only partial variables (reducing memory used)
fd = matfile(fd);

if strcmp(ga,'all')                 % To plot all grains
    ga = 1:size(thSB,1);
end

p = cell(length(ga),1);             % Vector with plots generated
c = 0;                              % Grain counter
for i = ga
    c = c + 1;
    gF = fd.gF(i,1);                % Partial loading of figures from fft decomposition
    gF = gF{1};
    gI = fp.gI(i,1);                % Partial loading of figures from preprocessing
    gI = gI{1};
 
    if any(strcmp(dispP,'gr'))                  % PLOT: Reduced greyscale image
        p{c}.gr = figure;
        imshow(gI.r(:,:,2))
    end
    
    if any(strcmp(dispP,'skel'))                % PLOT: Skeletonized grain
        p{c}.skel = figure;
        imshow(gF.If)
    end
    
    if any(strcmp(dispP,'thsearch'))            % PLOT: Theta mean contrast to find bright bands in FFT
        p{c}.thsearch = figure;
        dths = 180/length(smeanthi{i});
        ths = (0:dths:(180-dths));
        plot(ths,smeanthi{i},'LineWidth',1.5)
        hold on
        scatter(thF(i,:),pthF(i,:),600,'.','r')
        xlim([0,180])
        xticks(0:30:180)
        ylim([floor(min(smeanthi{i})-1),ceil(max(smeanthi{i})+1)])
        set(gca,'fontsize',16,'fontname','arial')
        xlabel('$\theta_s$ [$^{\circ}$]','interpreter','latex','fontsize',20)
        ylabel('$\overline{z_s}$','interpreter','latex','fontsize',20)
        box on
    end
    
    if any(strcmp(dispP,'fftbands'))            % PLOT: Bright bands in FFT
        p{c}.fftbands = figure;
        subplot(2,1,1)
        imshow(gF.Ifftd)
        subplot(2,1,2)
        imshow(gF.Ifftd)
        hold on
        [y0,x0] = deal(ceil((size(gI.g(:,:,1),1)+0.1)/2),ceil((size(gI.g(:,:,1),2)+0.1)/2));
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        xpi = rp.*cosd(thF(i,:)') + x0;
        ypi = rp.*sind(thF(i,:)') + y0;
        plot(xpi',ypi','LineWidth',1.5,'Color','c')
    end
    
    if any(strcmp(dispP,'cropped'))             % PLOT: Cropped FFTs
        p{c}.cropped = cell(size(gF.Ifftdc,3),1);
        for j = 1:size(gF.Ifftdc,3)
            p{c}.cropped{j} = figure;
            imshow(gF.Ifftdc(:,:,j))
        end
    end
    
    if any(strcmp(dispP,'icropped'))            % PLOT: Iffts from cropped FFTs
        p{c}.icropped = cell(size(gF.Iifftc,3),1);
        for j = 1:size(gF.Iifftc,3)
            p{c}.icropped{j} = figure;
            imshow(gF.Iifftc(:,:,j))
        end
    end
    
    if any(strcmp(dispP,'sbpeaks'))             % PLOT: Intensity spectra with slip band peaks identified
        p{c}.sbpeaks = figure;
        [y0,x0] = deal(ceil((size(gI.g(:,:,1),1)+0.1)/2),ceil((size(gI.g(:,:,1),2)+0.1)/2));
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        xp = rp.*cosd(thI(i,:)') + x0;
        yp = rp.*sind(thI(i,:)') + y0;
        for j = 1:sum(~isnan(thF(i,:)))
            subplot(sum(~isnan(thF(i,:))),2,2*(j-1)+1)
            imshow(gF.Iifftc(:,:,j))
            hold on
            plot(xp(j,:),yp(j,:),'LineWidth',1,'Color','w')
            scatter(pX{i}{j},pY{i}{j},15,'r','filled')
            subplot(sum(~isnan(thF(i,:))),2,2*(j-1)+2)
            plot((xp(j,:)-xp(j,1))/(xp(j,2)-xp(j,1))+1,S{i}(j,:))
            hold on
            scatter((pX{i}{j}-xp(j,1))/(xp(j,2)-xp(j,1))+1,sbI{i}{j},10,'r','filled')
        end
        sgtitle(['sbps = [',num2str(round(lambdameanp(i,~isnan(lambdameanp(i,:))),1)),'] pixel/SB'])
    end
    
    if any(strcmp(dispP,'sbI'))                 % PLOT: Slip bands identified in original image
        p{c}.sbI = figure;
        [y0,x0] = deal(ceil((size(gI.g(:,:,1),1)+0.1)/2),ceil((size(gI.g(:,:,1),2)+0.1)/2));
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        xp = rp.*cosd(thI(i,:)') + x0;
        yp = rp.*sind(thI(i,:)') + y0;
        imshow(gI.g(:,:,2))
        hold on
        for j = 1:sum(~isnan(thF(i,:)))
            plot(xp(j,:),yp(j,:),'LineWidth',1.5,'Color','c')
            scatter(pX{i}{j},pY{i}{j},25,'r','filled')
        end
    end
    
    if any(strcmp(dispP,'sbIadj'))                 % PLOT: Slip bands identified in original image
        p{c}.sbI = figure;
        [y0,x0] = deal(ceil((size(gI.g(:,:,1),1)+0.1)/2),ceil((size(gI.g(:,:,1),2)+0.1)/2));
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        xp = rp.*cosd(thI(i,:)') + x0;
        yp = rp.*sind(thI(i,:)') + y0;
        imshow(sqrt(gI.g(:,:,2)))
        hold on
        for j = 1:sum(~isnan(thF(i,:)))
            plot(xp(j,:),yp(j,:),'LineWidth',1.5,'Color','c')
            scatter(pX{i}{j},pY{i}{j},25,'r','filled')
        end
    end
    
    if any(strcmp(dispP,'sbIl'))                % PLOT: Slip bands identified in original image + SB at highest peak
        p{c}.sbIl = figure;
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        [X,Y] = meshgrid(1:size(gI.g(:,:,1),2),1:size(gI.g(:,:,1),1));
        imshow(gI.g(:,:,2))
        hold on
        for j = 1:sum(~isnan(thF(i,:)))
            if ~isempty(sbI{i}{j})
                for k = 1:length(sbI{i}{j})
                    xj = rp.*cosd(thSB(i,j)) + pX{i}{j}(k);
                    yj = rp.*sind(thSB(i,j)) + pY{i}{j}(k);
                    plot(xj,yj,'LineWidth',1.5,'Color','r')
                end
            end
        end
        scatter(X(gI.out),Y(gI.out),1,[1,1,1]*0,'.')
    end
    
    if any(strcmp(dispP,'sbprof'))              % PLOT: Comparison of slip bands profiles in orignal and IFFT images
        bpc = 8;                                    % Bands per row
        [y0,x0] = deal(ceil((size(gI.g(:,:,1),1)+0.1)/2),ceil((size(gI.g(:,:,1),2)+0.1)/2));
        rp = -(ceil(norm(size(gI.g(:,:,1))))/2):(ceil(norm(size(gI.g(:,:,1))))/2);
        xp = rp.*cosd(thI(i,:)') + x0;
        yp = rp.*sind(thI(i,:)') + y0;
        p{c}.sbprofG = cell(length(sbI{i}),1);
        p{c}.sbprofP = cell(length(sbI{i}),1);
        for j = 1:length(sbI{i})
            if ~isempty(sbI{i}{j})
                bj = length(sbI{i}{j});
                row = min(bpc,bj);
                col = floor(bj/bpc)+1;
                p{c}.sbprofG{j} = figure;
                imshow(gF.Iifftc(:,:,j))
                hold on
                plot(xp(j,:),yp(j,:),'LineWidth',1.5,'Color','c')
                scatter(pX{i}{j},pY{i}{j},25,'r','filled')
                
                p{c}.sbprofP{j} = figure;
                xnotnan = find(sum(~(isnan(Ssb0{i}{j}) | Ssb0{i}{j}==0),1));
                xlimj = xnotnan([1,end]);
                for k = 1:bj
                    subplot(row,col,k)
                    plot(1:length(Ssb0{i}{j}(1,:)),Ssb0{i}{j}(k,:))
                    hold on
                    plot(1:length(Ssb0{i}{j}(1,:)),Ssbifft{i}{j}(k,:))
                    xlim(xlimj)
                    xticks(diff(xlimj)*[0,1/4,1/2,3/4,1]+xlimj(1))
                    ylim([0,1])
                    yticks([0,0.5,1])
                    if k == 1
                        xticklabels({'0',[],[],[],[num2str(diff(xlimj)),' pDIC']})
                        xtickangle(0)
                        yticklabels({'0',[],'1'})
                        set(gca,'fontsize',9,'fontname','arial')
                    else
                        xticklabels([])
                        yticklabels([])
                    end
                    set(gcf,'color','w')
                    grid on
                    box on
                end
            end
        end
    end
end


end


