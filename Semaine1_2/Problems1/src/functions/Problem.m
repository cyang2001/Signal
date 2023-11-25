function Problem()
  [s1, Fs1] = audioread("MarteauPiqueur01.mp3");
  [s2, Fs2] = audioread("Jardin01.mp3");
  [s3, Fs3] = audioread("Jardin02.mp3");
  [s4, Fs4] = audioread("Ville01.mp3");
  Fs3
  signals = {s1, s2, s3, s4};
  Fs = {Fs1, Fs2, Fs3, Fs4};
  Ts = {1/Fs1, 1/Fs2, 1/Fs3, 1/Fs4};
  S = -48; % dBV
  G = 40; % dB
  P_SPL = 80; % dB SPL
  D_t = 1; % s
  windowSize = 37100; 
  sonStruct = struct('acceptable', struct(), 'penible', struct());
  for i = 1:length(signals)
    signal = signals{i};

    amplified_signal = FunctionAmplifier(signal, G);
    %figure;
    %numWindows = length(amplified_signal) - windowSize + 1; 
    %t = (0:numWindows-1) * Ts{i};  
    %plot(t,CalculateWindowedPowerdBm(amplified_signal, windowSize));
    %title(i);

    %peux etre changé
    stepSize = 37100; 
    dureeTemp1 = 0;  
    dureeTemp2 = 0;  
    Ts_i = Ts{i};
    sonTemp1 = [];
    sonTemp2 = []; 
    penibleStartTimes = [];
    penibleEndTimes = [];
    threshold_dBm = FunctionConvertSPLTodBM(P_SPL, S)+G;
    inHighSegment = false;
    segmentStartTime = 0;
    for n = 1:stepSize:length(amplified_signal) - windowSize + 1
        windowEnd = min(n + windowSize - 1, length(amplified_signal)); 
        currentWindow = amplified_signal(n:windowEnd);  
        current_dBm = CalculateWindowedPowerdBm(currentWindow, length(currentWindow));
        sonTemp2 = [sonTemp2 amplified_signal(n:windowEnd)]; 
        dureeTemp2 = dureeTemp2 + Ts_i*stepSize;
        if current_dBm >= threshold_dBm
            sonTemp1 = [sonTemp1 amplified_signal(n:windowEnd)]; 
            dureeTemp1 = dureeTemp1 + Ts_i*stepSize;
            inHighSegment = true;
            segmentStartTime = (n - 1 - stepSize/2) * Ts_i; 
        else
            if inHighSegment
                if dureeTemp1 >= D_t
                    segmentEndTime = (n - 1) * Ts_i; 
                    penibleStartTimes = [penibleStartTimes segmentStartTime(end)];
                    penibleEndTimes = [penibleEndTimes segmentEndTime];
                    sonStruct = FunctionSupport(i, dureeTemp1, sonTemp1, sonStruct, false, penibleStartTimes(end), penibleEndTimes(end));
                    sonTemp2 = sonTemp2(1:end-length(sonTemp1)+1); 
                    dureeTemp2 = dureeTemp2 - dureeTemp1; 
                    sonStruct = FunctionSupport(i, dureeTemp2, sonTemp2, sonStruct, true, 0, 0);
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
        sonStruct = FunctionSupport(i, dureeTemp1, sonTemp1, sonStruct, false, penibleStartTimes, penibleEndTimes);
    elseif dureeTemp2 > 0
        sonStruct = FunctionSupport(i, dureeTemp2, sonTemp2, sonStruct, true, 0, 0);
    end
end


sonStruct.acceptable
sonStruct.penible
for i = 1:length(signals)
    penibleSegments = [];
    for k = 1:length(sonStruct.penible.(['signal' num2str(i)]))
        penibleSegments = [penibleSegments, sonStruct.penible.(['signal' num2str(i)])(k).startTime, sonStruct.penible.(['signal' num2str(i)])(k).endTime];
    end
    
end

for i =1 : length(signals)
        penibleSegments = [];
        for k = 1:length(sonStruct.penible.(['signal' num2str(i)]))
            penibleSegments = [penibleSegments, sonStruct.penible.(['signal' num2str(i)])(k).startTime, sonStruct.penible.(['signal' num2str(i)])(k).endTime];
        end
        amplified_signal = FunctionAmplifier(signals{i}, G);
        t = (0:length(amplified_signal)-1) * Ts{i};  
        figure;
        subplot(2,1,1);
        plot(t,amplified_signal);
        hold on;
        for j = 1:2:length(penibleSegments)-1
            startTime = penibleSegments(j);
            endTime = penibleSegments(j+1);
            x = [startTime, endTime, endTime, startTime];
            y = [min(amplified_signal), min(amplified_signal), max(amplified_signal), max(amplified_signal)];
            patch(x, y, 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        numWindows = length(amplified_signal) - windowSize + 1; 
        subplot(2,1,2);
        t = (0:numWindows-1) * Ts{i};  
        plot(t,CalculateWindowedPowerdBm(signals{i}, windowSize));
        yline(FunctionConvertSPLTodBM(P_SPL,S), 'r--');
end
end

function sonStruct = FunctionSupport(i, dureeTemp, sonTemp, sonStruct, isAcceptable,penibleStartTimes, penibleEndTimes)
  signalN = sprintf('signal%d', i);
  duree = dureeTemp;
  son = sonTemp(1:end - 1);
  p_mW = FunctionCalculerPowerMeanmW(son);
  p_dBm = FunctionCalculerPowerMeandBM(son);
  rms = sqrt(p_mW);

  newSegment = struct('duree', duree, 'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms, 'son', son);
  newSegmentPenible = struct('duree', duree, 'puissance_mW', p_mW, 'puissance_dBm', p_dBm, 'Rms', rms, 'son', son ,'startTime', penibleStartTimes, 'endTime', penibleEndTimes);
  if isAcceptable
      if isfield(sonStruct.acceptable, signalN)
          sonStruct.acceptable.(signalN)(end+1) = newSegment;
      else
          sonStruct.acceptable.(signalN) = newSegment;
      end
  else
      if isfield(sonStruct.penible, signalN)
          sonStruct.penible.(signalN)(end+1) = newSegmentPenible;
      else
          sonStruct.penible.(signalN) = newSegmentPenible;
      end
  end
end


function signal_amplified = FunctionAmplifier(signal, G)
  signal_amplified = signal * 10 ^ (G / 20);
end

function P_dbm = FunctionConvertSPLTodBM(P_SPL, S)
  M_ref = 1;
  P_reference = 20 * 10 ^ (-6);
  P_dbm = 10 * log10((M_ref * 10 ^ (S / 20) * P_reference * 10 ^ (P_SPL / 20) * 10^2) ^ 2 * 1000);
end

function power_mean_dBm = FunctionCalculerPowerMeandBM(signal)
  power_mean_dBm = 10 * log10(FunctionCalculerPowerMeanmW(signal) / 0.001);
end

function power_mean_mW = FunctionCalculerPowerMeanmW(signal)
  power_mean_mW = mean(signal .^ 2) / 1000;
end

function p_dBm = CalculateWindowedPowerdBm(signal, windowSize)
    signalLength = length(signal);
    numWindows = signalLength - windowSize + 1;  
    p_mW = zeros(1, numWindows);
    
    for i = 1:numWindows
        windowStart = i;
        windowEnd = i + windowSize - 1;
        window = signal(windowStart:windowEnd);
        p_mW(i) = mean(window.^2);
    end

    p_dBm = 10 * log10(p_mW / 0.001);
end

