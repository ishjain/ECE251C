function [p0,theta] = win2Prototype(w,M,m)
if m == 1
    c = w(1:M/2); % The "cosine" terms
    s = flipud(w(M/2+1:M)); % The "sine" terms
    theta = atan2(s,c); % Estimate the angles based on the cosines and sines
    % Compute a new set of coefficients
    p0(1:M/2) = cos(theta);
    if isrow(p0)
        p0 = p0.';
    end
    p0(M/2+1:M) = flipud(sin(theta));
    p0(M+1:M+M/2) = flipud(p0(M/2+1:M));
    p0(M+M/2+1:2*M) = flipud(p0(1:M/2));
elseif m == 2
    for i = 1:M/2
        k = i-1;
        gk(i,:) = downsample(w,2*M,k);
        gMplusk(i,:) = downsample(w,2*M,M+k);
    end
    cd = gk(:,1)+gk(:,2);
    sd = gMplusk(:,2)+gMplusk(:,1);
    cs = gk(:,1)-gk(:,2);
    ss = gMplusk(:,2)-gMplusk(:,1);
    thetaDiff = atan2(sd,cd);
    thetaSum = atan2(ss,cs);
    theta0 = ((thetaSum+thetaDiff)/2);
    theta1 = ((thetaSum-thetaDiff)/2);
    theta = [theta0,theta1];
    % Convert the angles to polyphase components
    for i = 1:M/2
        gg(i,:) = [cos(theta0(i))*cos(theta1(i)), sin(theta0(i))*sin(theta1(i))];
        gg(2*M+1-i,:) = fliplr(gg(i,:));
        gg(M+i,:) = [-cos(theta0(i))*sin(theta1(i)), sin(theta0(i))*cos(theta1(i))];
        gg(M+1-i,:) = fliplr(gg(M+i,:));
    end
    % Assemble the polyphase components to form the prototype p0
    p0 = gg(:);
else
    error('m = %d is not supported!',m)
end