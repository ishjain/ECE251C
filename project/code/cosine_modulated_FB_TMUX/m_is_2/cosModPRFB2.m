% An FB from the Kaiser window prototype
clc
clf
close all
clearvars

%% Filter bank parameter definition
M = 16;
m = 2;
N = 2*m*M-1;

%% Prototype filter definition
p0 = [-0.00266329377333073;-0.00257376536834749;-0.00239790217528832;-0.00213924890860423;-0.00180191908350559;-0.00138862187506442;-0.000898205302745550;-0.000323359077284117;0.000350464425821993;0.00114407679066039;0.00208164164198208;0.00318697646723001;0.00448050335794765;0.00597762502128151;0.00768831219612399;0.00961707308660796;0.0117598471789568;0.0141090196831605;0.0166485345419061;0.0193479766640667;0.0221627884878250;0.0250328451787731;0.0278841640602675;0.0306345339098379;0.0332024522930040;0.0355170748055090;0.0375259916879871;0.0391983668997395;0.0405229489958718;0.0415024005862131;0.0421462459008232;0.0424644517025736;0.0424644517025736;0.0421462459008232;0.0415024005862131;0.0405229489958718;0.0391983668997395;0.0375259916879871;0.0355170748055090;0.0332024522930040;0.0306345339098379;0.0278841640602675;0.0250328451787731;0.0221627884878250;0.0193479766640667;0.0166485345419061;0.0141090196831605;0.0117598471789568;0.00961707308660796;0.00768831219612399;0.00597762502128151;0.00448050335794765;0.00318697646723001;0.00208164164198208;0.00114407679066039;0.000350464425821993;-0.000323359077284117;-0.000898205302745550;-0.00138862187506442;-0.00180191908350559;-0.00213924890860423;-0.00239790217528832;-0.00257376536834749;-0.00266329377333073].';

%% Filter bank definition
n = 0:N;
for idx = 1:M
    k = idx-1;
    thetak = (-1)^k*pi/4;
    h(idx,:) = 2*p0.*cos(pi/M*(k+1/2)*(n-N/2)+thetak);
    f(idx,:) = 2*p0.*cos(pi/M*(k+1/2)*(n-N/2)-thetak);
end

%% Test for PR
randseed = 0;
rng(randseed);
x = rand(1,100);
for idx = 1:M
    hout(idx,:) = conv(x,h(idx,:));
    DSout(idx,:) = downsample(hout(idx,:),M);
    USout(idx,:) = upsample(DSout(idx,:),M);
    fout(idx,:) = conv(USout(idx,:),f(idx,:));
end
xHat = sum(fout);
xHat = xHat*x(1)/xHat(64);