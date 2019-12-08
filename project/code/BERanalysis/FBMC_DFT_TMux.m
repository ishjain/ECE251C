



    W = exp(-1j*2*pi/N_SC);
    f0 = ones(N_SC+CP_LEN,1)/sqrt(N_SC); % Some arbitrary prototype filter

    S = ifft_in_mat;
    M=N_SC;
    FF = conj(dftmtx(M)); % The IDFT matrix
    repnum = ceil(length(f0)/M);
    Q = FF*S;
    % Q=ifft(S);
    Q = repmat(Q,repnum,1);
    Q = Q(end-N_SC-CP_LEN+1:end,:);
    Q = Q.*f0;
    move = N_SC+CP_LEN;
    for idx = 1:size(Q,2)
        i = idx-1;
        x4(idx,1:i*move+length(f0)) = [zeros(1,i*move),Q(:,idx).'];
    end
    tx_payload_vec = sum(x4);