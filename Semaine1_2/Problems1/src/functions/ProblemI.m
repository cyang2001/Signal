function ProblemI()
  [s1, Te1, s2, Te2, s3, Te3, s4, Te4] = FunctionF();
  signals = {s1, s2, s3, s4};
  Ts = {Te1, Te2, Te3, Te4};
  S = -48; % dBV
  G = 30; % dB
  P_SPL = 80; % dB SPL
  D_t = 1; % s

  for i = 1:length(signals)
      signal = signals{i};
      amplified_signal = FunctionAmplifier(signal, G);
      spl_signal = ConvertToSPL(amplified_signal, S);
      duree = length(signal); 
      duree = duree*Ts{i};
      power_mean_mW = mean(amplified_signal.^2)/1000;
      power_mean_dBm = 10 * log10(power_mean_mW / 0.001);
      rms = sqrt(power_mean_mW);

      if duree > D_t && max(spl_signal) >= P_SPL
          signal_class = 'son pénible';
      else
          signal_class = 'son acceptable';
      end

      resultats(i) = struct('duree', duree, 'power_mean_mW', power_mean_mW, ...
                            'power_mean_dBm', power_mean_dBm, 'rms', rms, ...
                            'signal_amplified', amplified_signal, 'signal_class', signal_class);
  end
  for i=1:length(resultats)
    fprintf('Signal %d:\n', i);
    fprintf('  Durée: %.4f s\n', resultats(i).duree);
    fprintf('  Puissance moyenne: %.4f mW (%.2f dBm)\n', resultats(i).power_mean_mW, resultats(i).power_mean_dBm);
    fprintf('  RMS: %.4f V\n', resultats(i).rms);
    fprintf('  Classe: %s\n', resultats(i).signal_class);
    fprintf('\n');
  end
end

function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10^(G/20);
end
function spl_signal = ConvertToSPL(signal, S)
  reference = 1; 
  spl_signal = 20 * log10(abs(signal) / reference) + (94 + S);
end
