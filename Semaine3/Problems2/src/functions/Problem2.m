function Problem2()
    % This is the path of the folder containing the audio files
    addpath ../../../Audios
    [s1, Fs1] = audioread('Pi_A_96K.wav');
    [s2, Fs2] = audioread('Pi_C_96K.wav');
    [s3, Fs3] = audioread('Vi_A3_96K.wav');
    [s4, Fs4] = audioread('Vi_C3_96K.wav');
    [s5, Fs5] = audioread('Vi_G4_96K.wav');
    [s6, Fs6] = audioread('Fl_A4_96K.wav');
    [s7, Fs7] = audioread('Fl_B3_96K.wav');
    % Calculate the average of the audio signal, to change it to mono channal
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

        % Detect note begin and end times
        [startTime, endTime] = DetectNoteTimes(p_W, Fs{i}, max(p_W));
        startTemp = [];
        endTemp = [];
        % Remove notes that are less than 1 second long
        for j = 1:length(startTime)
            if endTime(j) - startTime(j) >= 1
                startTemp = [startTemp, startTime(j)];
                endTemp = [endTemp, endTime(j)];
            end
        end

        disp(['The begining of this signal is :', num2str(startTemp), ' ending is :', num2str(endTemp)]);

        %Calculate average power  P
        avgPower = CalculateAveragePower(audioSignal{i}, Fs{i}, startTemp, endTemp);

        disp(['The average power of this note is :', num2str(avgPower), 'W']);

        %Calculate fundamental frequency
        f0 = AutocorrelationFundamentalFrequency(audioSignal{i}, Fs{i});

        disp(['The fundamental frequency of this note is :', num2str(f0), 'Hz']);

        %Calculate harmonics
        [fh, nh] = CalculateHarmonics(audioSignal{i}, Fs{i});
        disp(['The frequency of the highest harmonic of this note is :', num2str(fh), 'Hz']);
        disp(['The number of harmonics of this note is :', num2str(nh)]);

        subplot(2,1,2);
        N = length(audioSignal{i});
        Y = fft(audioSignal{i});
        P2 = abs(Y/N).^2; 
        f = (0:N-1)*(Fs{i}/N);
        plot(f(1:N/2+1), P2(1:N/2+1)); 
        xlabel('Frequency (Hz)');
        ylabel('Power/Frequency (dB/Hz)');
        title('Power Spectral Density of the Signal ');

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
  % Detect the start and end times of the notes in the given audio signal.

  % Threshold is 1% of the maximum power
  threshold = 0.01*max;

  % Find the indices where the signal is above the threshold
  aboveThreshold = audioSignal > threshold;

  % Make sure the signal is a column vector
  aboveThreshold = aboveThreshold(:);
  % Find the indices where the signal transitions from below to above the
  noteStartIndices = find(diff([0; aboveThreshold; 0]) == 1);
  noteEndIndices = find(diff([0; aboveThreshold; 0]) == -1) - 1;

  % Convert the indices to times
  startTime = noteStartIndices / fs;
  endTime = noteEndIndices / fs;
end

function avgPower = CalculateAveragePower(audioSignal, fs, startTime, endTime)
    startIndex = floor(startTime * fs);
    endIndex = ceil(endTime * fs);
    avgPower = mean(audioSignal(startIndex:endIndex).^2)/0.001;

end
function f0 = CalculateFundamentalFrequency(signal, fs)
  % Calculate the fundamental frequency of the given signal

    % Calculate the FFT of the signal
    N = length(signal);
    Y = fft(signal);
    % Calculate the bilateral spectrum P2
    P2 = abs(Y/N);
    % And then calculate the unilateral spectrum P1 based on P2 and the signal length N.
    P1 = P2(1:N/2+1);
    % Multiply the power by 2 to account for the negative frequencies
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(N/2))/N;
    % Find the frequency with the maximum amplitude
    [~, loc] = max(P1);
    f0 = f(loc);
end
function [fh, nh] = CalculateHarmonics(signal, fs)
  % Calculate the frequency of the highest harmonic and the number of harmonics
    N = length(signal);
    Y = fft(signal);
    P2 = abs(Y/N);
    P1 = P2(1:N/2+1);
    % Calculate the cumulative power and the total power
    totalPower = sum(P1.^2);
    cumulativePower = cumsum(P1.^2);
    % Find the frequency of the highest harmonic
    fhIndex = find(cumulativePower >= 0.9999 * totalPower, 1, 'first');
    f = fs*(0:(N/2))/N;
    fh = f(fhIndex);
    % Arbitrary threshold for harmonic counting
    nh = sum(P1 > 0.01 * max(P1));  
end

function f0_autocorr = AutocorrelationFundamentalFrequency(signal, fs)
  % Calculate the fundamental frequency of the given signal using the
  % autocorrelation method
    % Calculate the autocorrelation of the signal
    [autocorr, lags] = xcorr(signal, 'coeff');
    % Find the first positive peak of the autocorrelation
    positiveLags = lags >= 0;
    positiveAutocorr = autocorr(positiveLags);
    positiveLags = lags(positiveLags);
    [pks, locs] = findpeaks(positiveAutocorr, 'MinPeakProminence', 0.3, 'MinPeakDistance', fs/20);

    if ~isempty(locs)
        % Convert the lag to frequency
        peakLag = positiveLags(locs(1));
        f0_autocorr = fs / peakLag;
    else
        f0_autocorr = NaN;
    end
end

