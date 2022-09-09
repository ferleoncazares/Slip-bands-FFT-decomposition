%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to load a smoothed EBSD with specified 'grainId' property and
% remove the grain borders with a thickness of p pixels.
%
% Inputs:
% ebsd          Smoothed EBSD file with 'grainId' property
% pm             Number of pixels to remove from the grain borders 
%
% Outputs:
% gs            Surving grains [Id of surviving grain, pixels that survived]


function [g,gsbl,gAbl,gtop] = m_EBSD_gborderless(g,gs,gA,wm,Ablth,gtop)

g(gA < Ablth) = [];                         % Removing grains that do not pass the threshold
gs(gA < Ablth) = [];
gtop(gA < Ablth) = [];

gsbl = nan(size(g));                        % Surviving grains
gAbl  = nan(size(g));
for i = 1:length(gs)                        % Loop for each grain
    for j = 1:wm                            % Loop to eliminate pd pixels from the border
        k = boundary(g{i}.x,g{i}.y,1);      % Indices of points that correspond to the border
        g{i}(k) = [];                       % Discard the border
    end
    if length(g{i}) >= Ablth                % Grains above the threshold remaining pixels                     
        gsbl(i) = gs(i);                    % Index of grains ordered by size
        gAbl(i) = length(g{i});             % Surviving number of pixels
    end   
end

% Deleting grains (from all arrays) that did not survive the border reduction
g(isnan(gsbl(:)),:) = [];
gtop(isnan(gsbl(:)),:) = [];
gsbl(isnan(gsbl(:)),:) = [];
gAbl(isnan(gAbl(:)),:) = [];       
% ebsd_bl = cat(1,g{:});                    % Borderles EBSD

end

