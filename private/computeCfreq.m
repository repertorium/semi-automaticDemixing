function [X_ft,miditobins,muestrasmidi]=computeCfreq(x,fft_params,draw)
% [X_ft,miditobins]=computeCfreqSeq(x,fft_params,draw)
% Define the default fft_paramseters from the input signal
%
% Inputs:
% x  - signal
% fft_params - input parameters 
% draw - show process bar (default 0) 
%
% Outputs:
% X_ft - amplitude "MIDI" spectrogram (mnotesXframes)
% miditobins - Relation MIDI note & MIDI bins (2Xmnotes)
%
% Julio Carabias y Francisco Rodriguez. Fall 2012

if nargin<3,
    draw = 0;
end

% Inicializaciones
fs = fft_params.fs;
fftsize = fft_params.fftsize;
nframes = fft_params.nframes;
midi_inc = fft_params.midi_inc;
midi_min = fft_params.midi_min;
midi_max = fft_params.midi_max;

% Genero escala MIDI
miditobinskmin = zeros(1,(midi_max-midi_min+1)*midi_inc);
miditobinskmax = zeros(1,(midi_max-midi_min+1)*midi_inc);

for nota_midi=midi_min:midi_max,

    for midi_interval = 0:midi_inc-1,
        
        step_midi = 1/midi_inc;

        fmin=((2^(((nota_midi+midi_interval*step_midi-step_midi/2)-69)/12))*440);
        kmin=ceil(fmin/fs*(2*fftsize)+1);
        kmin=min(kmin,fftsize+1);

        fmax=((2^(((nota_midi+midi_interval*step_midi+step_midi/2)-69)/12))*440);
        kmax=fix(fmax/fs*(2*fftsize)+1);
        kmax=min(kmax,fftsize);

        miditobinskmin((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmin;
        miditobinskmax((nota_midi-midi_min)*midi_inc+midi_interval+1) = kmax;
        
    end;

end;

if midi_inc==1,
    miditobinskmax(miditobinskmax<miditobinskmin)=miditobinskmin(miditobinskmax<miditobinskmin);
else   
    indval = (miditobinskmax>=miditobinskmin);
    miditobinskmin = miditobinskmin(indval);
    miditobinskmax = miditobinskmax(indval);
    minkmin = min(miditobinskmin); % Lineal hasta el primer midi
    miditobinskmin = [1:minkmin-1 miditobinskmin];
    miditobinskmax = [1:minkmin-1 miditobinskmax];
end;

miditobinskmax(end) = fftsize+1;
muestrasmidi = length(miditobinskmin);
miditobins=[miditobinskmin;miditobinskmax];

% Inicializacion de variables de salida
X_ft=zeros(nframes,muestrasmidi);
X_ft = X_ft';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculo el spectrograma secuencialmente
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

windowsize = fft_params.windowsize;
hopsize = fft_params.hopsize;
nfft=2*fftsize;
win = sqrt(hanning(windowsize,'periodic'));  % Quitado SQRT a la ventana. FCO. 24/11/11
noverlap = windowsize-hopsize;

if ~any(any(imag(x)))    % x purely real
    if rem(nfft,2),    % nfft odd
        select = [1:(nfft+1)/2];
    else
        select = [1:nfft/2+1];
    end
else
    select = 1:nfft;
end

nx = length(x);
nwind = length(win);
if nx < nwind    % zero-pad x if it has length less than the win length
    x(end+1:nwind)=0;  nx=nwind;
end
x = x(:);     % make a column vector for ease later
win = win(:); % be consistent with data set

zp = zeros(nfft-nwind,1);
nwind2 = (nwind- mod(nwind,2))/2;
wzp = [win(nwind2+1:nwind);zp;win(1:nwind2)]; % !____!

%nframes = nframes;
nframes = fix((nx-noverlap)/(nwind-noverlap));
frame_index = 1 + (0:(nframes-1))*(nwind-noverlap);
if length(x)<(nwind+frame_index(nframes)-1)
    x(end+1:nwind+frame_index(nframes)-1) = 0;   % zero-pad x
end

% La puta draw
if draw,
    h = waitbar(0,'Please wait...','Name','Computing Spectrogram ...');
end;

% X = zeros(length(select),nframes);
for t=1:nframes,
    % Actualizo la draw
    if draw && mod(t,10)==0,
        waitbar(t/nframes,h,['Frame ' num2str(t) ' of ' num2str(nframes)]);
    end;
    
    xframe = x(frame_index(t):frame_index(t)+nwind-1); % extract frame of input data
    xzp = [xframe(nwind2+1:nwind);zp;xframe(1:nwind2)];
    xw = wzp .* xzp;
    Xframe = fft(xw); % FFT
    Xframe=Xframe(select);

    % Calculo del espectrograma en frecuencia midi
    for midi_index=1:muestrasmidi,
        
        kmin=miditobinskmin(midi_index);
        kmax=miditobinskmax(midi_index);
        
        if (kmax-kmin)==0,
            X_ft(midi_index,t) = abs(Xframe(kmin));
        elseif (kmax-kmin)>0,
            X_ft(midi_index,t) = sqrt(sum(abs(Xframe(kmin:kmax)).^2)); % / (kmax-kmin+1));
        end;
    end;
end;

% Cierro la puta draw
if draw
    close(h);
end;

return;