clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Credit:
% A detailed write-up of this example is available on the wiki:
% http://warpproject.org/trac/wiki/WARPLab/Examples/OFDM
%
% Copyright (c) 2015 Mango Communications - All Rights Reserved
% Distributed under the WARP License (http://warpproject.org/license)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all;

% Waveform params
N_OFDM_SYMS             = 100;         % Number of OFDM symbols
MOD_ORDER               = 4;           % Modulation order (2/4/16/64 = BSPK/QPSK/16-QAM/64-QAM)
TX_SCALE                = 1.0;          % Scale for Tx waveform ([0:1])

% OFDM params
N_SC                    = 64;                                     % Number of subcarriers/channels
CP_LEN                  = 16;                                     % Cyclic prefix length
INTERP_RATE             = 2;                                      % Interpolation rate (must be 2)
SAMP_FREQ           = 1e9;%20e6;


%% LTS for channel estimation

lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
SC_IND_PILOTS           = [8 22 44 58];                           % Pilot subcarrier indices
SC_IND_DATA             = [2:7 9:21 23:27 39:43 45:57 59:64];     % Data subcarrier indices
N_DATA_SYMS             = N_OFDM_SYMS * length(SC_IND_DATA);      % Number of data symbols (one per data-bearing subcarrier per OFDM symbol)
lts_t = ifft(lts_f, N_SC);

% Preamble contains 2.5 times lts for channel estimation
preamble = [ lts_t(N_SC/2+1:N_SC) lts_t lts_t];

%% Generate a payload of random integers
tx_data = randi(MOD_ORDER, 1, N_DATA_SYMS) - 1;

% Functions for data -> complex symbol mapping (like qammod, avoids comm toolbox requirement)
% These anonymous functions implement the modulation mapping from IEEE 802.11-2012 Section 18.3.5.8
modvec_bpsk   =  (1/sqrt(2))  .* [-1 1];
modvec_16qam  =  (1/sqrt(10)) .* [-3 -1 +3 +1];
modvec_64qam  =  (1/sqrt(43)) .* [-7 -5 -1 -3 +7 +5 +1 +3];

mod_fcn_bpsk  = @(x) complex(modvec_bpsk(1+x),0);
mod_fcn_qpsk  = @(x) complex(modvec_bpsk(1+bitshift(x, -1)), modvec_bpsk(1+mod(x, 2)));
mod_fcn_16qam = @(x) complex(modvec_16qam(1+bitshift(x, -2)), modvec_16qam(1+mod(x,4)));
mod_fcn_64qam = @(x) complex(modvec_64qam(1+bitshift(x, -3)), modvec_64qam(1+mod(x,8)));

% Map the data values on to complex symbols
switch MOD_ORDER
    case 2         % BPSK
        tx_syms = arrayfun(mod_fcn_bpsk, tx_data);
    case 4         % QPSK
        tx_syms = arrayfun(mod_fcn_qpsk, tx_data);
    case 16        % 16-QAM
        tx_syms = arrayfun(mod_fcn_16qam, tx_data);
    case 64        % 64-QAM
        tx_syms = arrayfun(mod_fcn_64qam, tx_data);
    otherwise
        fprintf('Invalid MOD_ORDER (%d)!  Must be in [2, 4, 16, 64]\n', MOD_ORDER);
        return;
end

% Reshape the symbol vector to a matrix with one column per OFDM symbol
tx_syms_mat = reshape(tx_syms, length(SC_IND_DATA), N_OFDM_SYMS);

% Define the pilot tone values as BPSK symbols
pilots = [1 1 -1 1].';

% Repeat the pilots across all OFDM symbols
pilots_mat = repmat(pilots, 1, N_OFDM_SYMS);

%% IFFT

% Construct the IFFT input matrix
ifft_in_mat = zeros(N_SC, N_OFDM_SYMS);

% Insert the data and pilot values; other subcarriers will remain at 0
ifft_in_mat(SC_IND_DATA, :)   = tx_syms_mat;
ifft_in_mat(SC_IND_PILOTS, :) = pilots_mat;

doFBMC=0;
if(doFBMC)
    %     take ifft_in_mat 64x100 for 64 channels and 100 time samples
    %     and return tx_payload_mat (80x100 for OFDM)
    W = exp(-1j*2*pi/N_SC);
    f0 = [1,2,3,4,3,5,1]; % Some arbitrary prototype filter
    %%--A channel synthesizer in direct form
    % Loop over each channel
    for idx = 1:N_SC
        k = idx-1;
        f(idx,:) = f0.*W.^(-k*[0:length(f0)-1]); % The filter
        w(idx,:) = upsample(ifft_in_mat(idx,:),N_SC); % The symbols after upsampling
        p(idx,:) = conv(w(idx,:),f(idx,:)); % The filtered symbols
    end
    tx_payload_mat = sum(p); % Sum over all channels
    
    % Reshape to a vector
    tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));
    
else
    %Perform the IFFT
    tx_payload_mat = ifft(ifft_in_mat, N_SC, 1);
    
    % Insert the cyclic prefix
    if(CP_LEN > 0)
        tx_cp = tx_payload_mat((end-CP_LEN+1 : end), :);
        tx_payload_mat = [tx_cp; tx_payload_mat];
    end
    % Reshape to a vector
    tx_payload_vec = reshape(tx_payload_mat, 1, numel(tx_payload_mat));
end


% Construct the full time-domain OFDM waveform
tx_vec = [preamble tx_payload_vec];

tx_vec_air=tx_vec;

%% Plot Results

figure;
subplot(2,1,1);
plot(real(tx_vec_air), 'b');
axis([0 length(tx_vec_air) -TX_SCALE TX_SCALE])
grid on;
title('Tx Waveform (I)');

subplot(2,1,2);
plot(imag(tx_vec_air), 'r');
axis([0 length(tx_vec_air) -TX_SCALE TX_SCALE])
grid on;
title('Tx Waveform (Q)');




