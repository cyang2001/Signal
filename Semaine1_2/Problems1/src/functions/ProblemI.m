function ProblemI()
  [s1, Te1, s2, Te2, s3, Te3, s4, Te4] = FunctionF();
  signals = {s1, s2, s3, s4};
  Ts = {Te1, Te2, Te3, Te4};
  S = -48; % dBV
  G = 30; % dB
  P_SPL = 80; % dB SPL
  D_t = 1; % s
  sonStruct = struct('acceptable', struct(), 'penible', struct());
  for i = 1:length(signals)
      signal = signals{i};
      amplified_signal = FunctionAmplifier(signal, G);
      dureeTemp = 0;
      n = 0;
      j = 0;
      sonTemp = [];
      for tempSignal = amplified_signal
          isLargerThandBm = false;
          n = n + 1;
          j = j + 1;
          sonTemp = [sonTemp, tempSignal];
          dureeTemp = dureeTemp + n*Ts{j};
          if FunctionCalculerPowerMeandBM(tempSignal) >= FunctionConvertSPLTodBM(P_SPL, S)
              isLargerThandBm = true;
          end
          if FunctionCalculerPowerMeandBM(tempSignal) < FunctionConvertSPLTodBM(P_SPL, S)
              if isLargerThandBm && dureeTemp < D_t
                FunctionSupport(i, dureeTemp, sonTemp, sonStruct, true);
                dureeTemp = 0;
                n = 0;
                sonTemp = [];
                sonTemp = [sonTemp, tempSignal];
              end
              if dureeTemp >= D_t && isLargerThandBm
                FunctionSupport(i, dureeTemp, sonTemp, sonStruct, false);
                dureeTemp = 0;
                n = 0;
                sonTemp = [];
                sonTemp = [sonTemp, tempSignal];
              end
          end
      end 

  end
  
end
function FunctionSupport(i, dureeTemp, sonTemp, sonStruct,isAcceptable)
  signalN = sprintf('signal%d', i);
  duree = dureeTemp;
  son = sonTemp(1:end-1);
  p_mW = FunctionCalculerPowerMeanmW(son);
  p_dBm = FunctionCalculerPowerMeandBM(son);
  rms = sqrt(p_mW);
  if isAcceptable
    sonStruct.acceptable.(signalN) = struct('duree',duree,'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms,'son',son);
  else
    sonStruct.penible.(signalN) = struct('duree',duree,'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms,'son',son);
  end
end
function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10^(G/20);
end
function P_dbm = FunctionConvertSPLTodBM(P_SPL, S)
  M_ref = 1;
  P_reference = 20*10^(-6);
  P_dbm = 10*log((M_ref*10^(S/20)*P_reference*10^(P_SPL/20)*30)^2*1000);
end
function power_mean_dBm = FunctionCalculerPowerMeandBM(signal)
  power_mean_dBm = 10 * log10(FunctionCalculerPowerMeanmW(signal) / 0.001);
end
function power_mean_mW = FunctionCalculerPowerMeanmW(signal)
  power_mean_mW = mean(signal.^2)/1000;
end