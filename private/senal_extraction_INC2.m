function [senales]=senal_extraction_INC2(x,param,Xe_ftj)
% [senales]=senal_extraction_INC(x,param,channel,A_ptj,S_fpj,M_js)
% Perform the signal source separation using Wiener Filtering.
%
% Further info can be found in:
% "Nonnegative signal factorization with learnt instrument models for sound 
%  source separation in close-microphone recordings"
%  JJ Carabias-Orti, M Cobos, P Vera-Candeas, FJ Rodr�guez-Serrano
%  EURASIP Journal on Advances in Signal Processing 2013 (1), 1-16
%
% Input arguments:
%	x --> input signal
%   params --> FFT parameters
%   channel --> Selected channel
%   A_ptj --> NMF Gains
%   S_fpj --> NMF Basis Functions
%   M_js --> Panning Matrix
%
% Output:
%   senales --> separated signals
%
% Francisco Rodriguez. January 2014
%

% Inicializaciones
j_max = size(Xe_ftj,3);

% Definicion de la salida
senales = zeros(length(x),j_max);
% 
% Xe_ftj = zeros([size(S_fpj,1),size(A_ptj,2),j_max]);
% for jj=1:j_max
%     Y_ft = S_fpj(:,:,jj) * A_ptj(:,:,jj);
%     Xe_ftj(:,:,jj) = Y_ft * M_js(jj,channel);
% end

%%% La separacion por casilla trama-midi
%%% Se multiplica por el valor de la trama con respecto al total
%%% Xe_ftj_relativo = Xe_ftj ./ repmat(sum(Xe_ftj,3)+eps,[1,1,j_max]);
%Xe_ftj_relativo = Xe_ftj ./ repmat(sum(Xe_ftj,3)+eps,[1,1,j_max]);
%Xe_ftj_relativo(repmat(sum(Xe_ftj,3)==0,1,j_max))=1/j_max;
denominador = (sum(Xe_ftj.^2,3));
Xe_ftj_relativo = (Xe_ftj.^2) ./ repmat(denominador+realmin,[1,1,j_max]);
% Xe_ftj_relativo(repmat(denominador==0,1,j_max))=1/j_max;
Xe_ftj_relativo(repmat(denominador==0,1,j_max))=0;

% Valores de configuracion
ventana = sqrt(hanning(param.windowsize,'periodic'));
% ventana = hanning(param.windowsize,'periodic');

% Espectrogramas
X=sg(x(:,1),2*param.fftsize,param.fs,ventana,param.windowsize-param.hopsize);
%X=sg_slow(x(:,1),2*param.fftsize,param.fs,ventana,param.windowsize-param.hopsize,1);
X_r = zeros(size(X,1),size(X,2),j_max);

% Calculo del espectrograma en frecuencia midi INCREMENTADA a 1/4 midi
for inc_midi=1:length(param.flog2bins)
    
    kmin=param.flog2bins(1,inc_midi);
    kmax=param.flog2bins(2,inc_midi); 

    if (kmax-kmin)>=0
        for jj=1:j_max
            X_r(kmin:kmax,1:size(Xe_ftj,2),jj) = X(kmin:kmax,1:size(Xe_ftj,2)) .* repmat(Xe_ftj_relativo(inc_midi,:,jj),kmax-kmin+1,1);
        end
    end 
end

% Para cada instrumento
for jj=1:j_max
    % Espectrograma inverso
    senal_act = invspecgram2(X_r(:,:,jj),2*param.fftsize,param.fs,ventana(:),param.windowsize-param.hopsize); % Ventana vector columna
    senales(1:length(senal_act),jj) = senal_act;
end

return;