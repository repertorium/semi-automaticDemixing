function res = calc_win_factor(M, WINDOW, NOVERLAP)
% res = calc_win_factor(NFRAMES, WINDOW, NOVERLAP)
%
% Compute window compensation factor for invspecgram2.
%

WINDOW = WINDOW(:);
WL = length(WINDOW);
SPACING = WL-NOVERLAP;

Tpad = NOVERLAP + M*SPACING;

Wind = zeros(Tpad,1);
WINDOW = WINDOW.^2;
for t = 1:M
    idx = ((t-1)*SPACING + (1:WL))';
    Wind(idx) = Wind(idx) + WINDOW;
end
Wind(1:round(WL/2)) = max(1,Wind(1:round(WL/2)));
Wind(end-round(WL/2):end) = max(1,Wind(end-round(WL/2):end));
res = max(Wind);
