function compute_Xfts_train(inputfolder)
% compute_Xfts_train - Multichannel spectrogram to learn panning.
%   The spectrogram (and other variables) is stored in file 
%   rX_fts_train.mat. This spectrogram can be used to learn the
%   panning matrix M_js.
%
% Syntax:   compute_Xfts_train(inputfolder)
%
% Inputs:
%    inputfolder - Folder contaning WAVs (44,1 kHz) and marker.txt
%
% Author: Pablo Cabanas
% email: pcabanas@ujaen.es
% Jan 2023


%% Max excerpt length for training (sec)
MAX_EXCERPT_LEN = 30;

%% FFT parameters
fs = 44100;
fft_params.hop = 0.01;       % Hop size in secons
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

%% Compute transform (rX_fts) with annotations for each interval
rX_fts = [];
intervals = [];
fprintf('Computing training spectrogram ....   0 %%');

for i = 1:length(markers)
    
    if markers(i).j == 0, continue; end
    
    rX_fts_chunk = [];
    
    inisample = markers(i).sample;
    endsample = Inf;
    if i < length(markers), endsample = markers(i+1).sample-1; end
    
    for s = 1:length(Dwav)
        [ry, fsori] = audioread(fullfile(inputfolder, Dwav(s).name), [inisample endsample]);
        if fsori ~= fs
            error('Sample-rate must be %d.', fs);
        end
        
        % Max excerpt lenght
        ry = ry(1:min(length(ry), MAX_EXCERPT_LEN*fs));
        
        % Freq transform
        rparam = fftParams(ry, fs, fft_params);
        rparam.midi_inc = 4;
        rX_fts_chunk(:,:,s) = computeCfreq(ry, rparam, 0); %#ok<AGROW>
    end
    
    % Annotation
    intervals(end+1).symbol = markers(i).symbol; %#ok<AGROW>
    intervals(end).iniframe = size(rX_fts, 2) + 1;
    intervals(end).endframe = size(rX_fts, 2) + size(rX_fts_chunk, 2);
    intervals(end).j        = markers(i).j;
    
    % Concatenate
    rX_fts = cat(2, rX_fts, rX_fts_chunk);
    
    fprintf('\b\b\b\b\b%3d %%', round(100*i/length(markers)) );
end
clear rX_fts_chunk ry j;
fprintf('\b\b\b\b\b%3d %%\n', 100);

%% NMF parameters
NMFparams.B             = 1.5;
NMFparams.NUM_MAX_ITER  = 50;
NMFparams.ALPHA         = 11.5129;
NMFparams.m_max         = 20;
NMFparams.j_max         = length(intervals);
NMFparams = NMF_setParams(rX_fts, rparam, NMFparams);
NMFparams.f_max         = size(rX_fts, 1);
NMFparams.t_max         = size(rX_fts, 2);
NMFparams.s_max         = size(rX_fts, 3);
NMFparams.p_max         = NMFparams.p_max * rparam.midi_inc;
NMFparams.k_max         = NMFparams.j_max * NMFparams.p_max;

save(fullfile(inputfolder, 'rX_fts_train.mat'), ...
    'rX_fts', 'intervals', 'NMFparams', 'Dwav', '-v7.3');
fprintf('File rX_fts_train.mat created\n');

return;
