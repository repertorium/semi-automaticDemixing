function [fft_params]=fftParams(x,fs,fft_params)
% [fft_params]=fft_params(x,fs,fft_params)
% Define the default fft_parameters from the input signal
%
% Inputs:
% x the sinal
% fs frequency sampling rate
%
% Outputs:
% fft_params - structure with default fft_paramss

fft_params.fs=fs;              
fft_params.ts=1/fft_params.fs;      

% Window and hop size in frames
% ReMAS
fft_params.hopsize = 570;%round(fft_params.hop*fs/2)*2;   % power of 2
fft_params.windowsize = 5700;%round(fft_params.window*fs/2)*2; % power of 2
% NO_ReMAS
% fft_params.hopsize = round(fft_params.hop*fs/2)*2;   % power of 2
% fft_params.windowsize = round(fft_params.window*fs/2)*2; % power of 2

% Update parameters (from power of 2)
fft_params.hop = fft_params.hopsize/fs;
fft_params.window = fft_params.windowsize/fs;

% Compute number of frames
fft_params.nframes = fix((length(x)-fft_params.windowsize)/fft_params.hopsize)+1;

% Log range (from midi low to high notes)
fft_params.midi_min=24;             
fft_params.midi_max=ceil(12*log2((fs/2)/440)+69);

% Numero de MIDI bins por cada nota
fft_params.midi_inc=1;

return;