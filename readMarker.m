function markers = readMarker(filename)
%readMarker - Read instruments marker file.
%
% Syntax:  markers = readMarker(filename)
%
% Inputs:
%    filename - Marker file (full path)
%
% Outputs:
%    markers - Struct with markers
%
%
% Authors: P. Cabanas-Molero (pcabanas@ujaen.es)
%          A.J. Munoz-Montoro (munozantonio@uniovi.es)
%          J.J. Carabias-Orti (carabias@ujaen.es)
% Last revision: Jan 2023


%% List of allowed instruments (2-character code)
instruments = {'PF', ...    % Piano
    'VN', ...               % Violin
    'VL', ...               % Viola
    'VC', ...               % Cello
    'CB', ...               % Contrabass
    'TR', ...               % Trumpet
    'HR', ...               % Horn
    'OB', ...               % Oboe
    'CL', ...               % Clarinet
    'FL', ...               % Flute
    'FG', ...               % Basson
    'TB', ...               % Trombone
    'TU', ...               % Tuba
    'AS', ...               % Alto Sax
    'TI', ...               % Timpani
    'TM'};                  % Tambour

%% Read file into memory
txt = fileread(filename);

%% Get markers
pat = '(?<hh>\d\d):(?<mm>\d\d):(?<ss>\d\d):(?<ms>\d\d\d)\s+(?<symbol>\w+)';
markers = regexp(txt, pat, 'names');

%% Determine time in samples and instrument index
fs = 44100;
j = 0;
for i = 1:length(markers)
    
    % Time of this marker in samples
    markers(i).sample = round( fs * ( ...
        str2double(markers(i).hh)*60*60 + ...
        str2double(markers(i).mm)*60 + ...
        str2double(markers(i).ss) + ...
        str2double(markers(i).ms)*0.001 ) ) + 1;
    
    % Assign an unique index j to the instrument in this marker
    if ~any(strcmp(instruments, markers(i).symbol(1:2)))
        markers(i).j = 0;   % not an instrument
        continue;
    end
    idx = find(strcmp({markers(1:i-1).symbol}, markers(i).symbol), 1, 'first');
    if isempty(idx)
        j = j+1;
        markers(i).j = j;   % new instrument
    else
        %markers(i).j = markers(idx).j;  % repeated instrument
        markers(i).j = 0;  % repeated instrument (ignore)
    end
end

return;