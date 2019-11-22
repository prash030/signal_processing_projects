% AlphaRatioCalculator.m (c) Prasanth "Prash" Ganesan
% Author: Prash Ganesan <prasganesan.pg@gmail.com>

function [Alpha_closed_open_R] = AlphaRatioCalculator(hdr_open,hdr_closed,EEG_open,EEG_closed)
  %AlphaRatioCalculator Calculates the alpha ratio of EEG signals

  Fs = hdr_open.frequency(1); % Fs is same for closed.
  L=length(EEG_open);

  %design an FIR highpass 
  Fstop = .5;
  Fpass = 1;
  Apass = 0.01;
  Astop = 80;
  param = fdesign.highpass(Fstop,Fpass,Astop,Apass,Fs);
  filt_FIR = design(param,'equiripple','SystemObject',true);
  EEG_open_hp_filt = filter(filt_FIR.Numerator,1,EEG_open);
  EEG_closed_hp_filt = filter(filt_FIR.Numerator,1,EEG_closed);

  %design an FIR lowpass 
  Fstop = 80;
  Fpass = 45;
  Apass = 0.01;
  Astop = 80;
  param = fdesign.lowpass('fp,fst,ap,ast',Fpass,Fstop,Apass,Astop,Fs);
  filt_FIR = design(param,'equiripple','SystemObject',true);
  EEG_open_lp_filt = filter(filt_FIR.Numerator,1,EEG_open_hp_filt);
  EEG_closed_lp_filt = filter(filt_FIR.Numerator,1,EEG_closed_hp_filt);

  FFT_EEG_open = abs(fft(EEG_open_lp_filt));
  FFT_EEG_closed = abs(fft(EEG_closed_lp_filt));

  %compare alpha energy
  s8 = L/Fs*8;
  s15 = L/Fs*15;
  alpha_open = sum(FFT_EEG_open(s8:s15))/sum(FFT_EEG_open)*100;
  alpha_closed = sum(FFT_EEG_closed(s8:s15))/sum(FFT_EEG_closed)*100;

  Alpha_closed_open_R = alpha_closed/alpha_open*100;
end

