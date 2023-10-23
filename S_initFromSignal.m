function S_initFromSignal(inputfolder, nbases)
% S_initFromSignal - Initialize instrument bases from signal
%   and saves them in inputfolder/S_pfj.mat.
%
% Syntax:   S_initFromSignal(inputfolder, nbases)
%
% Inputs:
%    inputfolder - Folder contaning WAVs, marker.txt, M_js.mat
%    nbases - Number of bases per instrument (200 by default)
%
% Author: Pablo Cabanas
% email: pcabanas@ujaen.es
% Jan 2023


%% FFT parameters
fs = 44100;
fft_params.hop = 0.01;       % Hop size in seconds
fft_params.window = 0.128;   % Window size
fft_params.fftsize = 8192;   % FFT size

%% Get WAV file names (audio channels)
Dwav = dir(fullfile(inputfolder, '*.wav'));
if isempty(Dwav)
    error('No wav files in this folder.'); 
end

%% Get markers
markerfile = fullfile(inputfolder, 'marker.txt');
if ~exist(markerfile, 'file')
    error('No marker file in this folder.');
end
markers = readMarker(markerfile);

%% Find beginning/ending of concert (tag 'ALL')
m = find(strcmpi('ALL', {markers.symbol}));
inisample = markers(m).sample;
endsample = Inf;
if m < length(markers)
    endsample = markers(m+1).sample - 1;
end

%% Load M_js
if ~exist(fullfile(inputfolder, 'M_js.mat'), 'file')
    main_M_train_manual(inputfolder);
end
load(fullfile(inputfolder, 'M_js.mat'), 'M_js', 'NMFparams');

%% Initialize S_pfj
if nargin < 2
    nbases = 200;
end
S_pfj = zeros(nbases, NMFparams.f_max, NMFparams.j_max);

% For each instrument, the model is the max value in each bin along
% its spot signal (during the concert interval)
fprintf('Initializing bases S_pfj ....   0 %%');
for jj = 1:NMFparams.j_max
    
    % Find spot mic
    [~, ss] = max(M_js(jj,:));
    
    % Read concert audio
    x = audioread(fullfile(inputfolder, Dwav(ss).name), ...
        [inisample endsample]);
    
    % Compute pattern
    param = fftParams(x, fs, fft_params);
    param.midi_inc = 4;
    X_ft = computeCfreq(x, param, 0);
    
    S_f = max(X_ft, [], 2);
    S_f = normalizeBeta(S_f, NMFparams.B);
    S_pfj(:,:,jj) = repmat(S_f', nbases, 1);
    
    fprintf('\b\b\b\b\b%3d %%', round(100*jj/NMFparams.j_max));
end
fprintf('\n');

%% Save S_pfj
save(fullfile(inputfolder, 'S_pfj.mat'), 'S_pfj');
fprintf('File S_pfj.mat created\n');

return;
