function ProblemI()
[s1, s2, s3, s4] = FunctionF();
signal = struct('signals', [s1, s2, s3, s4],'duree',[] ,'power_mean_mW', [], 'power_mean_dBm', [], 'rms', [], 'signal_amplified', [],'signal class', []);
S = -48; %S=-48dBV
G = 30; %G=30dB
P_SPL = 80; %P_SPL=30Db
D_t = 1; %D_t = 1s
signal_amplified = FunctionAmplifier(signal.signals, G);
signal.signal_amplified = signal_amplified;
for eachSignal = signal.signal_amplified
  length = length(eachSignal);
  signal.duree = [signal.duree, length];
  signal.power_mean_mW = [signal.power_mean_mW,mean(eachSignal.^2)];
  signal.power_mean_dBm = [signal.power_mean_dBm,10*log10(power_mean_mW)];
  signal.rms = [signal.rms,sqrt(power_mean_mW)];
  if length > 1 && 
end
end
function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10^(G/20);
end
