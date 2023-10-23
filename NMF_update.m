function [A_ptj, S_pfj, M_js] = NMF_update(v_cfsep, M_js, S_pfj, A_ptj, NMFparams)

%--------------------------------------------------------------------------
% Parameters
f_max = NMFparams.f_max;                % # of frequency points
p_max = NMFparams.p_max;                % # of notes
j_max = NMFparams.j_max;                % # of instruments
k_max = p_max * j_max;                  % # of bases

t_max = NMFparams.t_max;                % # of frames
s_max = NMFparams.s_max;                % # of channels

NUM_MAX_ITER = NMFparams.NUM_MAX_ITER;  % # of iterations
B = NMFparams.B;                        % Beta

%--------------------------------------------------------------------------
% Reconstruction
Y_fts = nmf_generateY(NMFparams, S_pfj, A_ptj, M_js);

%--------------------------------------------------------------------------
% Compute divergence between Y_fts (model) and v_cfsep (input)
D_B = zeros(NUM_MAX_ITER+1, 1);
if B == 1
    D_B(1) = sum(sum(sum((v_cfsep .* log((v_cfsep./(Y_fts+eps))+eps)) - v_cfsep + Y_fts)));
elseif B == 0
    D_B(1) = sum(sum(sum((v_cfsep./(Y_fts+eps)) - log((v_cfsep./(Y_fts+eps))+eps) - 1)));
else
    D_B(1) = sum(sum(sum((1/(B*(B-1))) * (v_cfsep.^B + (B-1)*(Y_fts+eps).^B- B*v_cfsep.*(Y_fts+eps).^(B-1)))));
end
fprintf(' Div (iter = %d)  =  %f\n', 0, D_B(1));

%--------------------------------------------------------------------------
% Iterative update
for ite = 1:NUM_MAX_ITER
    
    %----------------------------------------------------------------------
    % Update A_ptj
    
    % Combine M_js and S_pfj to compute W_fks (panned bases)
    S_fpj = permute(S_pfj, [2 1 3]);
    W_fpjs = zeros(f_max, p_max, j_max, s_max);
    for jj = 1:j_max
        for ss = 1:s_max
            W_fpjs(:,:,jj,ss) = S_fpj(:,:,jj) * M_js(jj,ss);
        end
    end
    W_fks = reshape(W_fpjs, f_max, k_max, s_max);
    clear S_fpj W_fpjs;
    
    Y_fts = Y_fts + eps;
    
    % Numerator gradient
    aux = (Y_fts .^ (B-2)) .* v_cfsep;
    KL_num = zeros(k_max, t_max);
    for ss = 1:s_max
        KL_num = KL_num + W_fks(:,:,ss)' * aux(:,:,ss); % Producto matricial para cada canal
    end
    
    % Denominator gradient
    aux = Y_fts .^ (B-1);
    KL_den = zeros(k_max, t_max);
    for ss = 1:s_max
        KL_den = KL_den + W_fks(:,:,ss)' * aux(:,:,ss); % Producto matricial para cada canal
    end
    clear aux;
    
    % Gradient
    delta = KL_num ./ (KL_den + eps);
    
    % Update A_ptj
    A_kt = reshape(permute(A_ptj, [1 3 2]), k_max, t_max);
    A_kt = A_kt .* delta;
    A_ptj = permute(reshape(A_kt, p_max, j_max, t_max), [1 3 2]);
    clear A_kt;
    
    % Reconstruction
    Y_fts = nmf_generateY(NMFparams, S_pfj, A_ptj, M_js);
    
    
    %----------------------------------------------------------------------
    % Compute divergence between Y_fts (model) and v_cfsep (input)
    if B==1
        D_B(ite+1) = sum(sum(sum((v_cfsep .* log((v_cfsep./(Y_fts+eps))+eps)) - v_cfsep + Y_fts)));
    elseif B==0
        D_B(ite+1) = sum(sum(sum((v_cfsep./(Y_fts+eps)) - log((v_cfsep./(Y_fts+eps))+eps) - 1)));
    else
        D_B(ite+1) = sum(sum(sum((1/(B*(B-1))) * (v_cfsep.^B + (B-1)*(Y_fts+eps).^B- B*v_cfsep.*(Y_fts+eps).^(B-1)))));
    end
    fprintf(' Div (iter = %d)  =  %f\n', ite, D_B(ite+1));
    
    
    
    %----------------------------------------------------------------------
    % Update S_pfj
    
    if isfield(NMFparams,'updateS') && NMFparams.updateS
        
        % Combine M_js and A_ptj to compute H_kts (panned gains)
        H_ptjs = zeros(p_max, t_max, j_max, s_max);
        for jj = 1:j_max
            for ss = 1:s_max
                H_ptjs(:,:,jj,ss) = A_ptj(:,:,jj) * M_js(jj,ss);
            end
        end
        H_kts = reshape(permute(H_ptjs, [1 3 2 4]), k_max, t_max, s_max);
        clear H_ptjs;
        
        Y_fts = Y_fts + eps;
        
        % Numerator gradient
        aux = (Y_fts .^ (B-2)) .* v_cfsep;
        KL_num = zeros(f_max, k_max);
        for ss = 1:s_max
            KL_num = KL_num + aux(:,:,ss) * H_kts(:,:,ss)'; % Producto matricial para cada canal
        end
        
        % Denominator gradient
        aux = Y_fts .^ (B-1);
        KL_den = zeros(f_max, k_max);
        for ss = 1:s_max
            KL_den = KL_den + aux(:,:,ss) * H_kts(:,:,ss)'; % Producto matricial para cada canal
        end
        clear aux;
        
        % Gradient
        delta = KL_num ./ (KL_den + eps);
        
        % Update S_pfj
        S_fk = reshape(permute(S_pfj, [2 1 3]), f_max, k_max);
        S_fk = S_fk .* delta;
        S_pfj = permute(reshape(S_fk, f_max, p_max, j_max), [2 1 3]);
        clear S_fk;
        
        % Reconstruction
        Y_fts = nmf_generateY(NMFparams, S_pfj, A_ptj, M_js);
        
    end


    %----------------------------------------------------------------------
    % Update M_js
    
    if (NMFparams.updateM)
        
        % Combine S_pfj and A_ptj to compute Y_ftj (spectrogram per instrument)
        Y_ftj = zeros(f_max, t_max, j_max);
        for jj = 1:j_max
            Y_ftj(:,:,jj) = S_pfj(:,:,jj)' * A_ptj(:,:,jj);
        end
        
        Y_fts = Y_fts + eps;
        
        % Numerator gradient
        aux = (Y_fts .^ (B-2)) .* v_cfsep;
        KL_num_M = zeros(j_max, s_max);
        for ff = 1:f_max
            Y = squeeze(Y_ftj(ff,:,:));
            Aux = squeeze(aux(ff,:,:));
            
            KL_num_M = KL_num_M + Y' * Aux; % Producto matricial para cada frecuencia
        end
        
        % Denominator gradient
        aux = Y_fts .^ (B-1);
        KL_den_M = zeros(j_max, s_max);
        for ff = 1:f_max
            Y = squeeze(Y_ftj(ff,:,:));
            Aux = squeeze(aux(ff,:,:));
            
            KL_den_M = KL_den_M + Y' * Aux; % Producto matricial para cada frecuencia
        end
        clear aux Y_ftj;
        
        % Update M_js
        delta = KL_num_M ./ (KL_den_M + eps);
        M_js = M_js .* delta;
        
        % Reconstruction
        Y_fts = nmf_generateY(NMFparams, S_pfj, A_ptj, M_js);
        
    end
    
    %----------------------------------------------------------------------
    % Compute divergence between Y_fts (model) and v_cfsep (input)
    if B==1
        D_B(ite+1) = sum(sum(sum((v_cfsep .* log((v_cfsep./(Y_fts+eps))+eps)) - v_cfsep + Y_fts)));
    elseif B==0
        D_B(ite+1) = sum(sum(sum((v_cfsep./(Y_fts+eps)) - log((v_cfsep./(Y_fts+eps))+eps) - 1)));
    else
        D_B(ite+1) = sum(sum(sum((1/(B*(B-1))) * (v_cfsep.^B + (B-1)*(Y_fts+eps).^B- B*v_cfsep.*(Y_fts+eps).^(B-1)))));
    end
    fprintf(' Div (iter = %d)  =  %f\n', ite, D_B(ite+1));

end
