function ProblemI()
[s1, s2, s3, s4] = FunctionF();
signal = [s1, s2, s3, s4];
S = -48; %S=-48dBV
G = 30; %G=30dB
P_SPL = 80; %P_SPL=30Db
D_t = 1; %D_t = 1s
signal_amplified = FunctionAmplifier(signal_amplified, G);
for eachSignal = signal_amplified
  length = length(eachSignal);
  power_mean_mW = mean(eachSignal.^2);
  power_mean_dBm = 10*log10(power_mean_mW);
  rms = sqrt(power_mean_mW);
end
end
function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10^(G/20);
end