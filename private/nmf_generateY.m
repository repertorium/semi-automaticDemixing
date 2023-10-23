function [Y_fts]=nmf_generateY(NMFparams,S_pfj,A_ptj,M_js)

% Parametros
f_max = 401;                % N_MIDI_SEP
j_max = NMFparams.j_max;    % SParam.N_INSTRUMENTS
t_max = NMFparams.t_max;    % Tamaño del bloque
s_max = NMFparams.s_max;    % Numero de canales

Y_fts=zeros(f_max,t_max,s_max); % Inicializo a ceros la matriz de 3 dimensiones

for jj = 1:j_max
    Y_ft = S_pfj(:,:,jj)' * A_ptj(:,:,jj); % Producto matricial para cada instrumento
    
    for ss = 1:s_max
        Y_fts(:,:,ss) = Y_fts(:,:,ss) + (Y_ft * M_js(jj,ss));        
    end
    
end
return;
