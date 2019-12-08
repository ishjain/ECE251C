clear
close all


x  = randi([0 15],5000,1);
y = qammod(x,16);

L = 4; %for upsampling and reconstruction
yPulse = rectpulse(y,L);

acpr = comm.ACPR(...
    'SampleRate', 3.84e6*8,...
    'MainChannelFrequency', 0,...
    'MainMeasurementBandwidth', 3.84e6,...
    'AdjacentChannelOffset', [-5e6 5e6],...
    'AdjacentMeasurementBandwidth', 3.84e6,...
    'MainChannelPowerOutputPort', true,...
    'AdjacentChannelPowerOutputPort', true);

figure;
subplot(211)
plot(abs(yPulse))
subplot(212)
plot(abs(fftshift(fft(yPulse))))

[ACPR,mainChnlPwr,adjChnlPwr] = acpr(yPulse)