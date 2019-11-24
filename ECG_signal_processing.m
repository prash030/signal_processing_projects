% Copyright (c) Prasanth "Prash" Ganesan
% Author email: <prasganesan.pg@gmail.com>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%---------------------------------------------------
clear
% Add the wfdb-app-toolbox-0-9-10 toolbox to the path
% MIT BIH NSR database
% The database has 18 long-term ECG recordings and are recorded from patients 
% from the Arrhythmia Laboratory at Boston's Beth Israel Hospital. 
% All patients had normal heart rhythm and the population include 5 men, 
% aged 26 to 45, and 13 women, aged 20 to 50.

% Import ecg
[t,ECG]=rdsamp('nsrdb/16265',[1],1000);
ann=rdann('nsrdb/16265','atr',[],1000);
fs=128;
time_req=3; % required time is 3 sec
samples = time_req*fs;
new_time = linspace(0,time_req,samples);
new_ecg=ECG(1:samples);
figure; plot(new_time,new_ecg);
title('3 sec of ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

% Plot annotations over that
hold on;
plot(ann(ann<=samples)./fs,new_ecg(ann(ann<=samples)),'r*');

%%--------------------------------------------
clear
time_req=2*3600; % 2 hours data
fs=128;
samples = time_req*fs;
[t,ECG]=rdsamp('nsrdb/16265',[1],samples);
ECG=ECG(:);
ann=rdann('nsrdb/16265','atr',[],samples);
ann=ann(:,1);
[HR_vec,tot_avg_HR,var_HR] = cal_HR(ECG,ann);

%%-------------------------------------------------------------------
clear
% MIT-BIH AFib database
% Atrial fibrillation is one of the life-threatening arrhythmia that is
% cahracterized by irregular and fast beating of the atria. It is caused by
% disorganization of the signals in the atria surface. The MIT-BIH afib
% database includes 25 long-term recordings of mostly paroxysmal afib
% patients. They are 10 hours ECG signal sampled at 250 samples/sec. The
% recordings were obtained from Boston's Beth Israel Hospital with the
% operating bandwidth of ECG monitors to be 0.1 to 40 Hz. The amplitude
% range of the recordings are -10 to +10 millivolts.

% Import and Plot signal
[t,ECG]=rdsamp('afdb/04015',[1],1000);
fs=250;
time_req=3; % required time is 3 sec
samples = time_req*fs;
new_time = linspace(0,time_req,samples);
new_ecg=ECG(1:samples);
figure; plot(new_time,new_ecg);
title('3 sec of Atrial Fibrillation ECG')
xlabel('Time (sec)')
ylabel('Amplitude')

% Disease description
% The disease present in the signal is atrial fibrillation. It can be
% clearly seen that the signals are morphologically inconsistent with no
% p-wave present. Although the R waves are not noisy, the time intervals 
% of R-waves are completely irregular and the
% baseline is fractionated.

%%---------------------------------------------------
clear
% Import and Plot signal
[t,ECG1]=rdsamp('ecgiddb/Person_01/rec_5',[1],4000);
figure; plot(t,ECG1);
title('Time domain of ECG-ID database signal')
xlabel('Time (sec)')
ylabel('Amplitude')

% Plot in Frequency domain
N_fft= length(ECG1);
L = length(ECG1);
fs = 500;
ECG1_fft = abs(fft(ECG1));
f = [0:(N_fft-1)]*fs/L;
figure; plot(f,abs(ECG1_fft))
xlabel('frequency (Hz)'); title('Frequqncy domain of ECG')

% Design the high pass filter
% From the frequency domain plot, we can see that the noise is present at
% around 1 to 2 Hz (see the high amplitude at the beginning). Se we need to
% pass from 10 Hz and stop till 5 Hz. The amplitude is also chosen till
% 150.
Fstop = 5;  % Stopband Frequency
Fpass = 10;    % Passband Frequency
Astop = 150;    % Stopband Attenuation (dB)
Apass = 0.01;      % Passband Ripple (dB)
Fs    = 500;    % Sampling Frequency

h = fdesign.highpass('fst,fp,ast,ap', Fstop, Fpass, Astop, Apass, Fs);
IIR_HPF = design(h, 'ellip', 'SystemObject', true);

% Design low pass filter
% From the frequency domain plot we can see that the noise is present at 50
% Hz. Hence we need to choose the stop band to include 50 Hz.
Fpass = 45;   % Passband Frequency
Fstop = 55;   % Stopband Frequency
Apass = 0.01;    % Passband Ripple (dB)
Astop = 150;  % Stopband Attenuation (dB)
Fs    = 500;  % Sampling Frequency

h = fdesign.lowpass('fp,fst,ap,ast', Fpass, Fstop, Apass, Astop, Fs);
IIR_LPF = design(h, 'ellip', 'SystemObject', true);

% View the filters
fvtool(IIR_HPF,'Fs',Fs); title('IIR HPF frequency response')
[b_HPF,a_HPF]=sos2tf(IIR_HPF.SOSMatrix,IIR_HPF.ScaleValues);

fvtool(IIR_LPF,'Fs',Fs); title('IIR LPF frequency response')
[b_LPF,a_LPF]=sos2tf(IIR_LPF.SOSMatrix,IIR_LPF.ScaleValues);

% Filtered signal
y_IIR_HPF = filter(b_HPF,a_HPF,ECG1);
figure; plot(t,ECG1); hold on; plot(t,y_IIR_HPF);
xlabel('time(sec)'); title('IIR HPF filtered response (Baseline wander has been removed)')
legend({'Original signal','Baseline wander Filtered Signal'});

y_IIR_HPF_LPF = filter(b_LPF,a_LPF,y_IIR_HPF);
figure; plot(t,ECG1); hold on; plot(t,y_IIR_HPF_LPF);
xlabel('time(sec)'); title('IIR HPF+LPF filtered response (Baseline wander + powerline noise both has been removed)')
legend({'Original signal','Baseline wander + Powerline noise Filtered Signal'});

% Frequency Domain of filtered signal
y_fft = abs(fft(y_IIR_HPF_LPF));
N_FFT = length(y_fft);
L = length(y_fft);
f = [0:(N_FFT-1)]*Fs/L;
figure; plot(f,y_fft);
xlabel('frequency (Hz)');
title('FFT of the designed Filter Response')

