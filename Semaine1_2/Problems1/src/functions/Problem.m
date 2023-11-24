function Problem()
  [s1, Te1, s2, Te2, s3, Te3, s4, Te4] = FunctionF();
  signals = {s1, s2, s3, s4};
  Ts = {Te1, Te2, Te3, Te4};
  S = -48; % dBV
  G = 40; % dB
  P_SPL = 80; % dB SPL
  D_t = 1; % s
  sonStruct = struct('acceptable', struct(), 'penible', struct());
  for i = 1:length(signals)
    signal = signals{i};
    amplified_signal = FunctionAmplifier(signal, G);
    dureeTemp1 = 0;  
    dureeTemp2 = 0;  
    Ts_i = Ts{i};
    sonTemp1 = [];
    sonTemp2 = []; 
    threshold_dBm = FunctionConvertSPLTodBM(P_SPL, S);
    inHighSegment = false;

    for n = 1:length(amplified_signal)
        current_dBm = FunctionCalculerPowerMeandBM(amplified_signal(n));
        sonTemp2 = [sonTemp2 amplified_signal(n)]; 
        dureeTemp2 = dureeTemp2 + Ts_i;
        if abs(current_dBm) >= threshold_dBm
            sonTemp1 = [sonTemp1 amplified_signal(n)]; 
            dureeTemp1 = dureeTemp1 + Ts_i;
            inHighSegment = true;
        else
            if inHighSegment
                if dureeTemp1 >= D_t
                    sonStruct = FunctionSupport(i, dureeTemp1, sonTemp1, sonStruct, false);
                    sonTemp2 = sonTemp2(1:end-length(sonTemp1)+1); 
                    dureeTemp2 = dureeTemp2 - dureeTemp1; 
                    if size(sonTemp2) == 0
                        sonStruct = FunctionSupport(i, dureeTemp2, sonTemp2, sonStruct, true);
                    end
                    sonTemp2 = [];
                    dureeTemp2 = 0;
                    dureeTemp2 = dureeTemp2 + Ts_i;
                end
                sonTemp1 = [];
                dureeTemp1 = 0;
                inHighSegment = false;
            end
        end
    end

    if dureeTemp1 >= D_t
        sonStruct = FunctionSupport(i, dureeTemp1, sonTemp1, sonStruct, false);
    elseif dureeTemp2 > 0
        sonStruct = FunctionSupport(i, dureeTemp2, sonTemp2, sonStruct, true);
    end
end



sonStruct.acceptable
sonStruct.penible
sonStruct.acceptable.signal1.duree
sonStruct.penible.signal1
end

function sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, isAcceptable)
  signalN = sprintf('signal%d', i);
  duree = dureeTemp;
  son = sonTemp(1:end - 1);
  p_mW = FunctionCalculerPowerMeanmW(son);
  p_dBm = FunctionCalculerPowerMeandBM(son);
  rms = sqrt(p_mW);

  newSegment = struct('duree', duree, 'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms, 'son', son);

  if isAcceptable
      if isfield(sonStruct.acceptable, signalN)
          sonStruct.acceptable.(signalN)(end+1) = newSegment;
      else
          sonStruct.acceptable.(signalN) = newSegment;
      end
      %fprintf("Le signal %d est acceptable\n", i);
  else
      if isfield(sonStruct.penible, signalN)
          sonStruct.penible.(signalN)(end+1) = newSegment;
      else
          sonStruct.penible.(signalN) = newSegment;
      end
      %fprintf("Le signal %d est penible\n", i);
  end
end


function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10 ^ (G / 20);
end

function P_dbm = FunctionConvertSPLTodBM(P_SPL, S)
  M_ref = 1;
  P_reference = 20 * 10 ^ (-6);
  P_dbm = 10 * log((M_ref * 10 ^ (S / 20) * P_reference * 10 ^ (P_SPL / 20) * 10^2) ^ 2 * 1000);
end

function power_mean_dBm = FunctionCalculerPowerMeandBM(signal)
  power_mean_dBm = 10 * log10(FunctionCalculerPowerMeanmW(signal) / 0.001);
end

function power_mean_mW = FunctionCalculerPowerMeanmW(signal)
  power_mean_mW = mean(signal .^ 2) / 1000;
end


