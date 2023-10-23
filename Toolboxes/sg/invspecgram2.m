function res = invspecgram2(B, NFFT, Fs, WINDOW, NOVERLAP)
% res = invspecgram2(B,NFFT,Fs,WINDOW,NOVERLAP)
% All parameters are like specgram...
% Written by Niels Henrik Pontoppidan
% Modified by pcabanas to work frame-by-frame (low memory usage)
% Use calc_win_factor to compensate window factor.
%

[N,M] = size(B);
Next = NFFT-N;

WINDOW = WINDOW(:);
WL = length(WINDOW);
SPACING = WL-NOVERLAP;

% Reserve output vector
Tpad = NOVERLAP + M*SPACING;
res = zeros(Tpad,1);

% Index of beginnings of frames
frames_index = 1 + (0:(M-1))*SPACING;

% Frame-by-frame
for t = 1:M
    % Compute IFFT
    Bfull = zeros(NFFT,1);
    Bfull(1:N) = B(:,t);
    Bfull(N+1:end) = conj(Bfull(Next+1:-1:2));
    idft = real(ifft(Bfull));
    idft = idft(1:WL);
    
    % Overlap-add
    if t==1
        res(frames_index(1):frames_index(1)+WL-1) = idft.*WINDOW;
    else
        res(frames_index(t):frames_index(t)+WL-1) = idft.*WINDOW + ...
            res(frames_index(t):frames_index(t)+WL-1);
    end
end

% Compensate window
% Wind = zeros(Tpad,1);
% WINDOW = WINDOW.^2;
% for t = 1:M
%     idx = ((t-1)*SPACING + (1:WL))';
%     Wind(idx) = Wind(idx) + WINDOW;
% end
% Wind(1:round(WL/2)) = max(1,Wind(1:round(WL/2)));
% Wind(end-round(WL/2):end) = max(1,Wind(end-round(WL/2):end));
% res = res./max(Wind);

