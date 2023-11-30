function Problem2()
    addpath ../../../Audios
    [s1, Fs1] = audioread('Pi_A_96K.wav');
    [s2, Fs2] = audioread('Pi_C_96K.wav');
    [s3, Fs3] = audioread('Vi_A3_96K.wav');
    [s4, Fs4] = audioread('Vi_C3_96K.wav');
    [s5, Fs5] = audioread('Vi_G4_96K.wav');
    [s6, Fs6] = audioread('Fl_A4_96K.wav');
    [s7, Fs7] = audioread('Fl_B3_96K.wav');
    audioSignal = {mean(s1,2), mean(s2,2), mean(s3,2), mean(s4,2), mean(s5,2), mean(s6,2), mean(s7,2)};
    Fs = {Fs1, Fs2, Fs3, Fs4, Fs5, Fs6, Fs7};
    
    for i=1:length(audioSignal)
        windowSize = 20000; 
        p_W = CalculateWindowedPowerSliding(audioSignal{i}, windowSize);
        numWindows = length(audioSignal{i}) - windowSize + 1; 
        t = (0:numWindows - 1) * (1 / Fs{i}); 
        figure;
        subplot(2,1,1);
        plot(t, p_W);
        xlabel('Time (s)');
        ylabel('W');
        title('Audio Signal in Time Domain');
        [startTime, endTime] = DetectNoteTimes(p_W, Fs{i}, max(p_W));
        avgPower = CalculateAveragePower(audioSignal{i}, Fs{i}, startTime, endTime);
        subplot(2,1,2);
        N = length(audioSignal{i});
        Y = fft(audioSignal{i});
        f = (0:N-1)*(Fs{i}/N);
        plot(f, Y);
        f0 = CalculateFundamentalFrequency(audioSignal{i}, Fs{i});
        [fh, nh] = CalculateHarmonics(audioSignal{i}, Fs{i});
    end
    
end

function p_W = CalculateWindowedPowerSliding(signal, windowSize)
    % Calculate the power of the given signal in dBm using a sliding window
    % approach. This function smoothens the power calculation and converts 
    % it to dBm format.
    
    % Length of the signal
    signalLength = length(signal);

    % Number of windows based on the signal length and window size
    numWindows = signalLength - windowSize + 1;  

    % Initialize an array to store the power values for each window
    p_mW = zeros(1, numWindows);

    % Calculate the power for the first window
    window = signal(1:windowSize);
    currentPower = mean(window.^2);  % Mean power of the first window
    p_mW(1) = currentPower;

    % Iteratively calculate the power for the remaining windows
    for i = 2:numWindows
        % Update the power by adding the new element and subtracting the old
        currentPower = currentPower - signal(i - 1)^2 / windowSize + signal(i + windowSize - 1)^2 / windowSize;
        p_mW(i) = currentPower;
    end
    p_W = p_mW/0.001;
    
end

function [startTime, endTime] = DetectNoteTimes(audioSignal, fs, max)
  threshold = 0.01*max;
  aboveThreshold = audioSignal > threshold;
  aboveThreshold = aboveThreshold(:);
  noteStartIndices = find(diff([0; aboveThreshold; 0]) == 1);
  noteEndIndices = find(diff([0; aboveThreshold; 0]) == -1) - 1;
  startTime = noteStartIndices / fs;
  endTime = noteEndIndices / fs;
end
function avgPower = CalculateAveragePower(audioSignal, fs, startTime, endTime)
    startIndex = floor(startTime * fs);
    endIndex = ceil(endTime * fs);
    avgPower = mean(audioSignal(startIndex:endIndex).^2)/0.001;

end
function f0 = CalculateFundamentalFrequency(signal, fs)
    N = length(signal);
    Y = fft(signal);
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(N/2))/N;
    [~, loc] = max(P1);
    f0 = f(loc);
end
function [fh, nh] = CalculateHarmonics(signal, fs)
    N = length(signal);
    Y = fft(signal);
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    totalPower = sum(P1.^2);
    cumulativePower = cumsum(P1.^2);
    fhIndex = find(cumulativePower >= 0.9999 * totalPower, 1, 'first');
    f = fs*(0:(N/2))/N;
    fh = f(fhIndex);
    nh = sum(P1 > 0.01 * max(P1));  % Arbitrary threshold for harmonic counting
end

