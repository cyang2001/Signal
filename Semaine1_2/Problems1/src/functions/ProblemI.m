function ProblemI()
    [s1, Te1, s2, Te2, s3, Te3, s4, Te4] = FunctionF();
    signals = {s1, s2, s3, s4};
    Ts = {Te1, Te2, Te3, Te4};
    S = -48; % dBV
    G = 40; % dB
    P_SPL = 80; % dB SPL
    D_t = 1; % s
    i = 1;
    sonStruct = struct('acceptable', struct(), 'penible', struct());
    signal = signals{i};
    amplified_signal = FunctionAmplifier(signal, G);
    dureeTemp = 0;
    n = 0;
    Ts_i = Ts{i};
    sonTemp = [];
    fprintf('Type of amplified_signal: %s\n', class(amplified_signal));
    fprintf('Size of amplified_signal: %s\n', mat2str(size(amplified_signal)));
    threshold_dBm = FunctionConvertSPLTodBM(P_SPL, S);
    isLargerThandBm = false;
    windowSize = 1/Te1; 
    stepSize = 10; 

    for n = 1:stepSize:length(amplified_signal) - windowSize + 1
      currentWindow = amplified_signal(n:n + windowSize - 1);
      current_dBm = FunctionCalculerPowerMeandBM(currentWindow); 
  
      if current_dBm >= threshold_dBm
          if ~isLargerThandBm
              isLargerThandBm = true;
              dureeTemp = 0;
              sonTemp = [];
          end
          sonTemp = [sonTemp, currentWindow];
          dureeTemp = dureeTemp + Ts_i * windowSize;
      else
          if isLargerThandBm
              if dureeTemp >= D_t
                  sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, false);
              else
                  sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, true);
              end
              isLargerThandBm = false;
              sonTemp = [];
              dureeTemp = 0;
          end
      end
  end
  if isLargerThandBm
      if dureeTemp >= D_t
          sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, false);
      else
          sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, true);
      end
  end
  

    sonStruct.acceptable
    sonStruct.penible
end

function sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, isAcceptable)
    signalN = sprintf('signal%d', i);
    duree = dureeTemp;
    son = sonTemp(1:end - 1);
    p_mW = FunctionCalculerPowerMeanmW(son);
    p_dBm = FunctionCalculerPowerMeandBM(son);
    rms = sqrt(p_mW);

    if isAcceptable
        sonStruct.acceptable.(signalN) = struct('duree', duree, 'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms, 'son', son);
        fprintf("Le signal %d est acceptable\n", i);
    else
        sonStruct.penible.(signalN) = struct('duree', duree, 'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms, 'son', son);
        fprintf("Le signal %d est penible\n", i);
    end

end

function signal_amplified = FunctionAmplifier(signal, G)
    signal_amplified = signal * 10 ^ (G / 20);
end

function P_dbm = FunctionConvertSPLTodBM(P_SPL, S)
    M_ref = 1;
    P_reference = 20 * 10 ^ (-6);
    P_dbm = 10 * log((M_ref * 10 ^ (S / 20) * P_reference * 10 ^ (P_SPL / 20) * 30) ^ 2 * 1000);
end

function power_mean_dBm = FunctionCalculerPowerMeandBM(signal)
    power_mean_dBm = 10 * log10(FunctionCalculerPowerMeanmW(signal) / 0.001);
    fprintf('Power mean dBm: %d\n', power_mean_dBm);
end

function power_mean_mW = FunctionCalculerPowerMeanmW(signal)
    power_mean_mW = mean(signal .^ 2) / 1000;
    fprintf('Power mean mW: %d\n', power_mean_mW);
end
