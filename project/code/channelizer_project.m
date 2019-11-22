offsets = [-40,-30,-20,10,15,25,35,-15];
sinewave = dsp.SineWave('ComplexOutput',true,'Frequency',...
    offsets+(-375:125:500),'SamplesPerFrame',800);

channelizer = dsp.Channelizer('StopbandAttenuation',140);
synthesizer = dsp.ChannelSynthesizer('StopbandAttenuation',140);
spectrumAnalyzer =  dsp.SpectrumAnalyzer('ShowLegend',true,'NumInputPorts',...
    2,'ChannelNames',{'Input','Output'},'Title','Input and Output of QMF');

for i = 1:5000
    x = sum(sinewave(),2);
    y = channelizer(x);
    v = synthesizer(y);
    spectrumAnalyzer(x,v)
end