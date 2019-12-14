function [factor,val] = diffFromPowComp(p0,M,m)
% Extract the kth and the (M+k)th polyphase components out of the 2M
% polyphase components
for i = 1:M/2
    k = i-1;
    gk(i,:) = downsample(p0,2*M,k);
    gMplusk(i,:) = downsample(p0,2*M,M+k);
end
if m == 1
    % Normalize the coefficients
    factor = sqrt(1/(gk(1)^2+gMplusk(1)^2));
    gk = gk*factor;
    gMplusk = gMplusk*factor;

    % Compute the expressions
    for i = 1:M/2
        val(i) = abs(1-gk(i)^2-gMplusk(i)^2);
    end
elseif m == 2
    % Normalize the coefficients
    factor = sqrt(1/((gk(1,1)+gk(1,2))^2+(gMplusk(1,1)+gMplusk(1,2))^2));
    gk = gk*factor;
    gMplusk = gMplusk*factor;

    % Compute the expressions
    for i = 1:M/2
        val(i,1) = abs(1-(gk(i,1)+gk(i,2))^2-(gMplusk(i,1)+gMplusk(i,2))^2);
        val(i,2) = abs(1-(gk(i,1)-gk(i,2))^2-(-gMplusk(1,1)+gMplusk(1,2))^2);
    end
else
    error('m = %d is not supported!',m)
end