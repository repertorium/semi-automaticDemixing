function main_train(initfolder)
% MAIN_TRAIN - Multichannel source separation model initialization.
%
%   MAIN_TRAIN(INITFOLDER) initializes the panning matrix M_js and
%   the bases S_pfj of the model. The directory INITFOLDER must contain
%   the WAV files along with an annotation file called marker.txt. This
%   file must indicate: 1) the instruments in the mixture and 2) temporal
%   segments where each instrument plays solo. See the URMP folder for
%   an example marker.txt file.
%
%   Input:
%      INITFOLDER - folder containing WAVs and marker.txt
%
%   Example:
%      main_train('./URMP/31_Slavonic_tpt_tpt_hn_tbn/train');

%   Author: P. Cabanas-Molero (pcabanas@ujaen.es)
%   Last Revised: Jan 2023


rng('shuffle');
warning off;

%% Add to path
addpath(genpath('Toolboxes'));

%% Some params
initS = 'MODEL';    % SIGNAL or MODEL

%% Train M_js
M_train_manual(initfolder);

%% Initialize S_pfj
switch initS
    case 'SIGNAL'
        S_initFromSignal(initfolder);
    case 'MODEL'
        S_initFromModels(initfolder);
    otherwise
        S_initFromModels(initfolder);
end

end
