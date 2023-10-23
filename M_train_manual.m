function M_train_manual(inputfolder)
% M_train_manual - Initialize panning matrix and saves result
%   in inputfolder/M_js.mat.
%
% Syntax:   M_train_manual(inputfolder)
%
% Inputs:
%    inputfolder - Folder contaning WAVs (44,1 kHz) and marker.txt
%
% Author: Pablo Cabanas
% email: pcabanas@ujaen.es
% Jan 2023


%% Read training spectrogram
inputfile = fullfile(inputfolder, 'rX_fts_train.mat');
if ~exist(inputfile, 'file')
    compute_Xfts_train(inputfolder);
end
load(inputfile, 'rX_fts', 'intervals', 'NMFparams', 'Dwav');

%% Initialize M_js
M_js = ones(NMFparams.j_max, NMFparams.s_max);

%% Estimate panning coefficients for each instrument
for jj = 1:NMFparams.j_max
    % Compute energy in each channel
    iniframe = intervals(jj).iniframe;
    endframe = intervals(jj).endframe;
    En = squeeze( sum(rX_fts(:,iniframe:endframe,:).^2, [1 2]) )';
    
    % Panning coefficients
    M_js(jj,:) = sqrt(En / max(En));
end

%% Save M_js
save(fullfile(inputfolder, 'M_js.mat'), ...
    'M_js', 'intervals', 'NMFparams', 'Dwav', '-v6');
fprintf('File M_js.mat created\n');

return;
