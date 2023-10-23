function separate_blocks(inputfolder, initfolder, outputfolder, NMFparams)
% separate_blocks - Multichannel source separation based on
%   panning (block-by-block processing).
%
% Syntax:   separate_blocks(inputfolder, initfolder, outputfolder, NMFparams)
%
% Inputs:
%    inputfolder  - Folder contaning WAVs (44,1 kHz)
%    initfolder   - Folder containing M_js and S_pfj
%    outputfolder - Output folder
%    NMFparams    - Struct with parameters
%
% Authors: Pablo Cabanas (pcabanas@ujaen.es)
% Jan 2023


%% FFT parameters
fs = 44100;
fft_params.hop     = 0.01;      % Hop size in seconds
fft_params.window  = 0.128;     % Window size
fft_params.fftsize = 8192;      % FFT size

%% Load data
load('miditobins_sep.mat', 'miditobins_sep');

%% Get wav file names (audio channels)
Dwav = dir(fullfile(inputfolder, '*.wav'));
if isempty(Dwav)
    error('No wav files in this folder.');
end

%% Load M_js
load(fullfile(initfolder, 'M_js.mat'), 'M_js');

%% Load instrument bases S_pfj
load(fullfile(initfolder, 'S_pfj.mat'), 'S_pfj');

%% Analysis parameters
info     = audioinfo(fullfile(inputfolder, Dwav(1).name));
NSAMPLES = info.TotalSamples;
rparam   = fftParams(zeros(NSAMPLES,1), fs, fft_params);
NFRAMES  = rparam.nframes;
noverlap = rparam.windowsize - rparam.hopsize;

% Compute window compensantion factor
ventana   = sqrt(hanning(rparam.windowsize, 'periodic'));
winfactor = calc_win_factor(NFRAMES, ventana, noverlap);
clear ventana;

% Overlap segment
wavoverlp = zeros(noverlap, NMFparams.j_max, NMFparams.s_max);

%% NMF parameters
BLOCKSIZE       = min(NMFparams.BLOCKSIZE, NFRAMES);
NMFparams.t_max = BLOCKSIZE;
NMFparams.p_max = size(S_pfj, 1);

%% Create output WAV files
outwavs = createOutputFiles(outputfolder, {Dwav(:).name}, NMFparams.j_names);

% Write WAV header
for s = 1:NMFparams.s_max
    for j = 1:NMFparams.j_max
        writeWavHeader(outwavs{j,s}, NFRAMES*rparam.windowsize+noverlap, fs);
    end
end

%% Block-by-block processing
iniframe = 1;
endframe = BLOCKSIZE;

while iniframe <= NFRAMES
    inisample = (iniframe-1)*rparam.hopsize + 1;
    endsample = (endframe-1)*rparam.hopsize + rparam.windowsize;
    
    fprintf('Processing frames (%d-%d)/%d (%.2f %%) ...\n', ...
        iniframe, endframe, NFRAMES, endframe*100/NFRAMES);
    
    %% Read audio segment and compute rX_fts
    ry = [];
    rX_fts = [];
    for s = 1:NMFparams.s_max
        ry(:,s) = audioread(fullfile(inputfolder, Dwav(s).name), ...
            [inisample endsample]);                         %#ok<AGROW>
        
        rparam = fftParams(ry, fs, fft_params);             % Update # of frames
        rparam.midi_inc = 4;
        rX_fts(:,:,s) = computeCfreq(ry(:,s), rparam, 0);   %#ok<AGROW>
    end
    NMFparams.t_max = size(rX_fts, 2);
    
    % Normalize rX_fts
    rX_fts = normalizeMatrixBeta(rX_fts, NMFparams.B);
    
    %% Initialize A_ptj
    A_ptj = ones(NMFparams.p_max, NMFparams.t_max, NMFparams.j_max);
    
    %% Estimate A_kt and S_pfk
    [A_ptj, S_pfj_out] = NMF_update(rX_fts, M_js, S_pfj, A_ptj, NMFparams);
    
    %% Separate sources
    % Estimate spectrogram per instrument
    Y_ftj = zeros(NMFparams.f_max, NMFparams.t_max, NMFparams.j_max);
    for j = 1:NMFparams.j_max
        Y_ftj(:,:,j) = squeeze(S_pfj_out(:,:,j))' * A_ptj(:,:,j);
    end
    
    rparam.flog2bins = miditobins_sep;
    
    % Pan instruments to each channel and synthesize
    sources = zeros(endsample-inisample+1, NMFparams.j_max, NMFparams.s_max);
    for s = 1:NMFparams.s_max
        Y_ftj_pan = 0 * Y_ftj;
        for j = 1:NMFparams.j_max
            Y_ftj_pan(:,:,j) = squeeze(Y_ftj(:,:,j)) * M_js(j,s);
        end
        %Y_ftj_pan = postproc_mask(rX_fts, Y_ftj_pan, M_js, s);
        
        sources(:,:,s) = senal_extraction_INC2(ry(:,s), rparam, Y_ftj_pan);
    end
    
    % Overlap-add
    sources(1:noverlap,:,:) = sources(1:noverlap,:,:) + wavoverlp;
    wavoverlp(:,:,:)        = sources(end-noverlap+1:end,:,:);
    
    % Save output waveforms for this block
    for s = 1:NMFparams.s_max
        for j = 1:NMFparams.j_max
            samples = sources(1:end-noverlap, j, s) / winfactor;
            writeWavSamples(outwavs{j,s}, samples);
        end
    end
    
    iniframe = iniframe + BLOCKSIZE;
    endframe = endframe + BLOCKSIZE;
    endframe = min(endframe, NFRAMES);
end

return;


%--------------------------------------------------------------------------

function writeWavHeader(wavfile, nsamples, fs)

h = fopen(wavfile, 'w');

% Sizes (assume 16-bit samples)
DataSize = nsamples*2;
HeaderSize = 44;
ChunkSize = uint32(DataSize + HeaderSize - 8);
Subchunk2Size = uint32(DataSize + HeaderSize - 44);

% RIFF Chunk Descriptor
fwrite(h, 'RIFF', 'uchar');         % RIFF Header Magic header
fwrite(h, ChunkSize, 'uint32');     % RIFF Chunk Size
fwrite(h, 'WAVE', 'uchar');         % WAVE Header

% "fmt" sub-chunk
fwrite(h, 'fmt ', 'uchar');         % FMT header
fwrite(h, uint32(16), 'uint32');    % Size of the fmt chunk
fwrite(h, uint16(1), 'uint16');     % Audio format 1=PCM, 6=mulaw, 7=alaw,
                                    % 257=IBM Mu-Law, 258=IBM A-Law, 259=ADPCM
fwrite(h, uint16(1), 'uint16');     % Number of channels 1=Mono 2=Stereo
fwrite(h, uint32(fs), 'uint32');    % Sampling Frequency in Hz
fwrite(h, uint32(fs*2), 'uint32');  % Bytes per second
fwrite(h, uint16(2), 'uint16');     % 2=16-bit mono, 4=16-bit stereo
fwrite(h, uint16(16), 'uint16');    % Number of bits per sample

% "data" sub-chunk
fwrite(h, 'data', 'uchar');                 % "data"  string
fwrite(h, uint32(Subchunk2Size), 'uint32'); % Sampled data length

fclose(h);

return;

%--------------------------------------------------------------------------

function writeWavSamples(wavfile, samples)

h = fopen(wavfile, 'a');
samples(samples > 1)  =  1;
samples(samples <-1)  = -1;
samples = int16( (samples+1)*65535/2 - 65536/2 );
fwrite(h, samples, 'int16');
fclose(h);

return;

%--------------------------------------------------------------------------

function outwavs = createOutputFiles(outdir, wavnames, instnames)

mkdir(outdir);

s_max = length(wavnames);
j_max = length(instnames);

outwavs = cell(j_max, s_max);

% Create subdir for each channel
for s = 1:s_max
    [~, chname, ~] = fileparts(wavnames{s});
    mkdir(fullfile(outdir, chname));
    
    % Create file for each instrument
    for j = 1:j_max
        outwavs{j,s} = fullfile(outdir, chname, [instnames{j} '.wav']);
        h = fopen(outwavs{j,s}, 'w');
        fclose(h);
    end
end

return;

%--------------------------------------------------------------------------

function X = normalizeMatrixBeta(X, B)

norma = sum(X(:).^B);
norma = (norma.^(1/B)) / ((B*(B-1)).^(1/B));
X = X / norma;

return;
