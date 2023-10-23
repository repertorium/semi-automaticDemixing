function main_separate(inputfolder, initfolder, outputfolder)
% MAIN_SEPARATE - Multichannel source separation.
%
%   MAIN_SEPARATE(INFOLDER, INITFOLDER, OUTFOLDER) separates instrument
%   sources in a multichannel mixture using a panning-based NMF model.
%
%   Input:
%      INFOLDER - input folder with the mixture (WAV files)
%      INITFOLDER - folder with model parameters (M_js and S_pfj)
%      OUTFOLDER - output folder
%
%   Example:
%      infolder = './URMP/31_Slavonic_tpt_tpt_hn_tbn/mixtures';
%      initfolder = './URMP/31_Slavonic_tpt_tpt_hn_tbn/train';
%      outfolder = './OUTPUT';
%      main_separate(infolder, initfolder, outfolder);

%   Author: P. Cabanas-Molero (pcabanas@ujaen.es)
%   Last Revised: Jan 2023


rng('shuffle');
warning off;

%% Add to path
addpath(genpath('Toolboxes'));

%% USER DEFINED PATHS
outputfolder = fullfile(outputfolder, ['date_ ' datestr(now,30)]);

%% USER DEFINED OPTIONS
load(fullfile(initfolder, 'M_js.mat'), 'NMFparams', 'intervals');
NMFparams.updateM      = false;                 % Update M
NMFparams.updateS      = true;                  % Update S (bases)
NMFparams.BLOCKSIZE    = 2048;                  % Block size (frames)
NMFparams.NUM_MAX_ITER = 025;                   % # of iterations
NMFparams.j_names      = {intervals(:).symbol};

fprintf('Running source separation\n');
fprintf('... with BLOCKSIZE = %d frames\n', NMFparams.BLOCKSIZE);
fprintf('... with %d iterations\n', NMFparams.NUM_MAX_ITER);
fprintf('... with %d instruments and %d channels\n', ...
    NMFparams.j_max, NMFparams.s_max);
if NMFparams.updateS
    fprintf('... with bases update\n');
else
    fprintf('... with fixed bases\n');
end

%% Run separation
separate_blocks(inputfolder, initfolder, outputfolder, NMFparams);

end
