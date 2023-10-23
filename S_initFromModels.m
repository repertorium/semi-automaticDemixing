function S_initFromModels(inputfolder)
% S_initFromModels - Initialize instrument bases from models
%   and saves them in inputfolder/S_pfj.mat.
%
% Syntax:   S_initFromModels(inputfolder)
%
% Inputs:
%    inputfolder - Folder contaning WAVs (44,1 kHz) and marker.txt
%
% Author: Pablo Cabanas
% email: pcabanas@ujaen.es
% Jan 2023


ruta_bases = './INSTRUMENT_BASES/';

%% Load M_js file
if ~exist(fullfile(inputfolder, 'M_js.mat'), 'file')
    main_M_train_manual(inputfolder);
end
load(fullfile(inputfolder, 'M_js.mat'), 'intervals', 'NMFparams');

%% Load instrument bases S_pfj
S_pfj = zeros(NMFparams.p_max, NMFparams.f_max, NMFparams.j_max);
for jj = 1:NMFparams.j_max
    instcode = intervals(jj).symbol(1:2);
    load(fullfile(ruta_bases, [instcode '.mat']), 'S_pf');
    S_pfj(:,:,jj) = S_pf;
end

%% Save
save(fullfile(inputfolder, 'S_pfj.mat'), 'S_pfj', '-v6');
fprintf('File S_pfj.mat created\n');

return;
