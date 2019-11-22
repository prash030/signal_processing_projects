

%% Problem 1
clear; clc
[hdr_open,EEG_open] = edfread('S001R01.edf');
[hdr_closed,EEG_closed] = edfread('S001R02.edf');

%(a)
[Alpha_closed_open_R] = AlphaRatioCalculator(hdr_open,hdr_closed,EEG_open(62,:),EEG_closed(62,:));

%(b)
for channel=1:size(EEG_open,1)
    [Alpha_closed_open_R] = AlphaRatioCalculator(hdr_open,hdr_closed,EEG_open(channel,:),EEG_closed(channel,:));
    Alpha_closed_open_R_all(channel)=Alpha_closed_open_R;
end

%(c)
xlabs=hdr_open.label;
figure; plot([1:65],Alpha_closed_open_R_all); set(gca,'XTick',[1:2:65],'XTickLabel',{xlabs{1:2:65}})
xlabel('Channel names'); ylabel('Alpha ratio')
title('Alpha ratio of all channels')

%(d) 
% From the plot, we can see that the ratio is increasing for channels cp6
% to F5 and then decreasing for Ft7 to Tp7. The highest ratio is for channel Po8 
% and the lowest ratio is for channel T9
%% Problem 2
clear; clc; 
[hdr_closed,EEG_closed] = edfread('S001R02.edf');

samples = 1000;
t=[0:samples]./samples;
EEG_closed_Oz = EEG_closed(62,1:samples);

% (a)
[imf,ort,nbits] = emd(EEG_closed_Oz);

% (b)
N_fft= length(imf);
L = length(imf);
fs = 160;
f = [0:(N_fft-1)]*fs/L;
for temp=1:size(imf,1)
    figure; plot(f,abs(fft(imf(temp,:))))
    xlabel('frequency (Hz)');ylabel('Magnitude'); 
    title(['FFT of IMF ' num2str(temp)])
end

% (c)
% Since the frequecy range of the alpha band is 7.5 to 12.5 Hz, the IMFs
% 1 and 2 matches the frequency. SO IMF1 and 2 correspond to the alpha
% rhythm

% (d) and (e)
estimated = sum(imf(1:2,:));
figure; plot(estimated); xlabel('Samples'); ylabel('Sum of IMF')
title('Estimated alpha rhythm signal')

% (f)
figure; plot(EEG_closed_Oz - estimated); xlabel('Samples'); ylabel('Amplitude')
title('EEG without alpha rhythm')

%% Problem 3
clear;clc;
% Generate white gaussian noise
time_req = 10; % seconds
fs = 100; %Hz
samples = time_req*fs; % seconds
new_time = linspace(0,time_req,samples);
noise = randn(1,samples); % AWGN
b=1;
a=[1 -0.3 0.7];
filt_noise = filter(b,a,noise);

% (a) Noise and PSD
figure; subplot(2,2,1)
plot(new_time,noise); xlabel('samples'); ylabel('amplitude');
title('White Gaussian Noise')
subplot(2,2,2)
periodogram(noise); xlabel('Normalized Frequency'); ylabel('Power');
title('PSD of White Gaussian Noise')
subplot(2,2,3)
plot(new_time,filt_noise); xlabel('samples'); ylabel('amplitude');
title('Output of the system')
subplot(2,2,4)
periodogram(filt_noise);  xlabel('Normalized Frequency'); ylabel('Power');
title('PSD of the output')

% (b) Magnitude TF of the system
M = fvtool(b,a);

% (c) The basic operation of the filter is multiplication in the frequency
% domain and since the shape of the filter (see the magnitude response) is
% like a bell curve, the output PSD also takes the same shape of the
% filter.
% White gaussian noise basically means that it has all the frequencies
% present in it and hence the power will be almost zero with some
% fluctuations. The filter removed some components of the noise and because
% of that the power in the system goes above zero which is seen in the PSD.

% (d) AR 2 Model
[a_AR, b_AR] = aryule(filt_noise,2);
[H,W]=freqz(sqrt(b_AR),a_AR);
figure; periodogram(filt_noise); hold on
plot(W/pi,20*log10(2*abs(H)/(2*pi)),'r')
title('PSD and its AR2 model');
xlabel('Normalized Frequency'); ylabel('Power');
legend({'PSD of the output','AR2 Model'})

% (e) 
error_a = abs(a_AR - a)
error_b = abs(b - sqrt(b_AR))

% From the error vales we can see that the AR model is very close to the
% original one in terms of the parameters.

[H_org,W_org] = freqz(b,a);
figure; periodogram(filt_noise); hold on
plot(W/pi,20*log10(2*abs(H)/(2*pi)),'r')
plot(W_org/pi,20*log10(2*abs(H_org)/(2*pi)),'g')
title('Comparison of AR2 model and the original magnitude');
xlabel('Normalized Frequency'); ylabel('Power');
legend({'PSD of the output','AR2 Model','Original magnitude'})

% From the comparison figure, we can see that the red curve is close
% to the green one meaning that the AR2 model almost successfully modelled the
% response.

%% Problem 4
clear; clc;
rest = [1 0.83 0.47 0.08 -0.22 -0.37 -0.39 -0.26];
Fatigue = [1 0.9 0.66 0.36 0.07 -0.17 -0.32 -0.37];

% (a) Problem solved by hand. See the attatched page at the end.

% Rest
D=solve('1+(a_1*0.83)+(a_2*0.47)-(b_1)^2=0','0.83+(a_1)+(a_2*0.83)=0','0.47+(a_1*0.83)+(a_2)=0', 'a_1', 'a_2', 'b_1');
a_1_rest = unique(double(D.a_1))
a_2_rest = unique(double(D.a_2))
b_1_square_rest = unique(abs(double(D.b_1)))^2

% Fatigue
D=solve('1+(a_1*0.9)+(a_2*0.66)-(b_1)^2=0','0.9+(a_1)+(a_2*0.9)=0','0.66+(a_1*0.9)+(a_2)=0', 'a_1', 'a_2', 'b_1');
a_1_fatigue = unique(double(D.a_1))
a_2_fatigue = unique(double(D.a_2))
b_1_square_fatigue = unique(abs(double(D.b_1)))^2

% (b)
rest_fft = abs(fft(rest));
Fatigue_fft = abs(fft(Fatigue));

N_fft= length(rest_fft);
L = N_fft;
fs = 1;
f = [0:(N_fft-1)]*fs/L;
figure; subplot(2,1,1); plot(f,rest_fft)
xlabel('frequency (Hz)'); ylabel('Power'); title('PSD before fatigue')
subplot(2,1,2); plot(f,Fatigue_fft)
xlabel('frequency (Hz)'); ylabel('Power'); title('PSD during fatigue')

% (c)
% From the PSD, both before and after fatigue shows similar characteristics
% except for the maxima and the minima of the power that decreases for the
% fatigue case which makes sense when compared to the lower values of
% parameters in the AR2 model.

%% Problem 5
clear; clc;
load hrv1
figure; plot(hrv1); title('Heart rate data')
xlabel('Samples'); ylabel('Heart rate signal')

N=length(hrv1);
MDL=[];
for P=2:20
    AR_M = ar(hrv1,P);
    ep = AR_M.NoiseVariance;
    MDL(P-1) = (N*log10(ep))+(P*log10(N));
end

figure; plot([2:20],MDL); xlabel('P - model order'); ylabel('MDL');
title('Finding the optimal order of AR model for hrv1 signal')

% For this data the best AR model is of order 12 since it gives the minimum
% MDL score.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For EEG1 in problem 1
[~,EEG1] = edfread('S001R01.edf');
EEG1=EEG1(62,:);
N=length(EEG1);
MDL=[];
for P=2:30
    AR_M = ar(EEG1,P);
    ep = AR_M.NoiseVariance;
    MDL(P-1) = (N*log10(ep))+(P*log10(N));
end

figure; plot([2:30],MDL); xlabel('P - model order'); ylabel('MDL');
title('Finding the optimal order of AR model for EEG1')

% For this data the best AR model is of order 10 since it gives the minimum
% MDL score.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For EEG2 in problem 1
[~,EEG2] = edfread('S001R02.edf');
EEG2=EEG2(62,:);
N=length(EEG2);
MDL=[];
for P=2:30
    AR_M = ar(EEG2,P);
    ep = AR_M.NoiseVariance;
    MDL(P-1) = (N*log10(ep))+(P*log10(N));
end

figure; plot([2:30],MDL); xlabel('P - model order'); ylabel('MDL');
title('Finding the optimal order of AR model for EEG2')

% For this data the best AR model is of order 25 since it gives the minimum
% MDL score.
