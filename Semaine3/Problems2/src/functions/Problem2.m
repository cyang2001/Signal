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
    audioSignal = {mean(s1, 2), mean(s2, 2), mean(s3, 2), mean(s4, 2), mean(s5, 2), mean(s6, 2), mean(s7, 2)};
    Fs = {Fs1, Fs2, Fs3, Fs4, Fs5, Fs6, Fs7};
    audioName = {'Pi_A_96K', 'Pi_C_96K', 'Vi_A3_96K', 'Vi_C3_96K', 'Vi_G4_96K', 'Fl_A4_96K', 'Fl_B3_96K'};
    notesPiano = ["A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#"];
    notesViolin = ["G", "G#", "A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#"];
    notesFlute = ["C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#", "A", "Bb", "B"];
    octaves = ['0', '1', '2', '3', '4', '5', '6', '7'];
    [tablePiano, tableViolin, tableFlute] = GenerateFrequencyTables();
    
    for i = 1:length(audioSignal)
        windowSize = 20000;
        p_W = CalculateWindowedPowerSliding(audioSignal{i}, windowSize);
        numWindows = length(audioSignal{i}) - windowSize + 1;
        t = (0:numWindows - 1) * (1 / Fs{i});
        figure;
        subplot(2, 1, 1);
        plot(t, p_W);
        xlabel('Time (s)');
        ylabel('W');
        title('Audio Signal in Time Domain');
        hold on;
        % Detect note begin and end times
        [startTemp, endTemp] = DetectNoteTimes(p_W, Fs{i}, max(p_W));
        for k = 1:length(startTemp)
          halfWindowSize = floor(windowSize / 2);
          startIndex = max(1, (startTemp(k)*Fs{i} - halfWindowSize));
          endIndex = min(length(p_W), (endTemp(k)*Fs{i} + halfWindowSize));
          plot(t(startIndex:endIndex), p_W(startIndex:endIndex), 'r');
        end
        
        hold off;
        disp(['The begining of this signal is :', num2str(startTemp), ' ending is :', num2str(endTemp)]);

        % Calculate average power  P
        avgPower_dBm = CalculateAveragePower(audioSignal{i}, Fs{i}, startTemp, endTemp);

        disp(['The average power of this note is :', num2str(avgPower_dBm), 'dBm']);

        % Calculate fundamental frequency
        f0 = AutocorrelationFundamentalFrequency(audioSignal{i}, Fs{i});

        disp(['The fundamental frequency of this note is :', num2str(f0), 'Hz']);
        % Only for test
        if i == 1 || i == 2
            tolerance = 2;
            [octave, note] = ind2sub(size(tablePiano), find(abs(tablePiano - f0) <= tolerance));
            disp(['The note of this frequency is :', char(notesPiano(note)), octaves(octave)]);
        end

        if i == 3 || i == 4 || i == 5
            tolerance = 4;
            [octave, note] = ind2sub(size(tableViolin), find(abs(tableViolin - f0) <= tolerance));
            disp(['The note of this frequency is :', char(notesViolin(note)), octaves(octave)]);
        end

        if i == 6 || i == 7
            tolerance = 2;
            [octave, note] = ind2sub(size(tableFlute), find(abs(tableFlute - f0) <= tolerance));
            disp(['The note of this frequency is :', char(notesFlute(note)), octaves(octave)]);

        end

        %Calculate harmonics
        [fh, nh] = CalculateHarmonics(audioSignal{i}, Fs{i}, f0);
        disp(['The frequency of the highest harmonic of this note is :', num2str(fh), 'Hz']);
        disp(['The number of harmonics of this note is :', num2str(nh)]);

        subplot(2, 1, 2);
        N = length(audioSignal{i});
        Y = fft(audioSignal{i});
        P2 = abs(Y / N) .^ 2;
        f = (0:N - 1) * (Fs{i} / N);
        plot(f(1:floor(N / 2 + 1)), P2(1:floor(N / 2 + 1)));
        xlabel('Frequency (Hz)');
        ylabel('Power/Frequency (W/Hz)');
        title('Power Spectral Density of the Signal ');
        frame = getframe(gcf);
        im = frame2im(frame);

        imwrite(im, ['../../results/' , audioName{i} , '.png']);
        %instrument = ClassifyInstrument(CalculateADSR(p_W, Fs{i}, max(p_W)), CalculateSpectralCentroid(p_W, Fs{i}), CalculateZeroCrossingRate(p_W));
        %disp([instrument]);
    end

end


% CalculateWindowedPowerSliding calculates the windowed power of a signal using a sliding window approach.
% The function takes in a signal and a window size as input and returns the power values for each window.
% The power values are normalized by dividing them by 0.001.
%
% Inputs:
%   - signal: The input signal for which the windowed power needs to be calculated.
%   - windowSize: The size of the sliding window.
%
% Output:
%   - p_W: The power values for each window, normalized by dividing them by 0.001.
%
% Example usage:
%   signal = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
%   windowSize = 3;
%   p_W = CalculateWindowedPowerSliding(signal, windowSize);
function p_W = CalculateWindowedPowerSliding(signal, windowSize)

    % Length of the signal
    signalLength = length(signal);

    % Number of windows based on the signal length and window size
    numWindows = signalLength - windowSize + 1;

    % Initialize an array to store the power values for each window
    p_mW = zeros(1, numWindows);

    % Calculate the power for the first window
    window = signal(1:windowSize);
    currentPower = mean(window .^ 2); % Mean power of the first window
    p_mW(1) = currentPower;

    % Iteratively calculate the power for the remaining windows
    for i = 2:numWindows
        % Update the power by adding the new element and subtracting the old
        currentPower = currentPower - signal(i - 1) ^ 2 / windowSize + signal(i + windowSize - 1) ^ 2 / windowSize;
        p_mW(i) = currentPower;
    end

    p_W = p_mW / 0.001;

end

function [startTemp, endTemp] = DetectNoteTimes(audioSignal, fs, max)
    % Detect the start and end times of the notes in the given audio signal.

    % Threshold is 1% of the maximum power
    threshold = 0.01 * max;

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
    startTemp = [];
    endTemp = [];
    % Remove notes that are less than 1 second long
    for j = 1:length(startTime)

        if endTime(j) - startTime(j) >= 1
            startTemp = [startTemp, startTime(j)];
            endTemp = [endTemp, endTime(j)];
        end

    end

end


% CalculateAveragePower calculates the average power of an audio signal within a specified time range.
% 
% Inputs:
%   - audioSignal: the input audio signal
%   - fs: the sampling frequency of the audio signal
%   - startTime: the start time of the time range (in seconds)
%   - endTime: the end time of the time range (in seconds)
%
% Output:
%   - avgPower_dBm: the average power of the audio signal within the specified time range (in dBm)
%
% Example usage:
%   avgPower = CalculateAveragePower(audioSignal, fs, 0.5, 1.5);
function avgPower_dBm = CalculateAveragePower(audioSignal, fs, startTime, endTime)
  startIndex = floor(startTime * fs);
  endIndex = ceil(endTime * fs);
  avgPower_W = mean(audioSignal(startIndex:endIndex) .^ 2); 
  avgPower_dBm = 10 * log10(avgPower_W / 0.001); 
end



% CalculateFundamentalFrequency calculates the fundamental frequency of a given signal.
% 
% Inputs:
%   - signal: The input signal for which the fundamental frequency needs to be calculated.
%   - fs: The sampling frequency of the signal.
%
% Output:
%   - f0: The fundamental frequency of the signal.
%
% Example usage:
%   signal = sin(2*pi*100*(0:0.001:1));
%   fs = 1000;
%   f0 = CalculateFundamentalFrequency(signal, fs);
function f0 = CalculateFundamentalFrequency(signal, fs)
  % Calculate the fundamental frequency of the given signal

  % Calculate the FFT of the signal
  N = length(signal);
  Y = fft(signal);
  % Calculate the bilateral spectrum P2
  P2 = abs(Y / N);
  % And then calculate the unilateral spectrum P1 based on P2 and the signal length N.
  P1 = P2(1:N / 2 + 1);
  % Multiply the power by 2 to account for the negative frequencies
  P1(2:end - 1) = 2 * P1(2:end - 1);
  f = fs * (0:(N / 2)) / N;
  % Find the frequency with the maximum amplitude
  [~, loc] = max(P1);
  f0 = f(loc);
end


% CalculateHarmonics calculates the harmonics of a given signal.
% 
% Inputs:
%   - signal: the input signal
%   - fs: the sampling frequency of the signal
%   - f0: the fundamental frequency of the signal
%
% Outputs:
%   - fh: the frequencies of the detected harmonics
%   - nh: the number of detected harmonics
%
% Example usage:
%   signal = sin(2*pi*1000*(0:1/fs:1));
%   fs = 10000;
%   f0 = 1000;
%   [fh, nh] = CalculateHarmonics(signal, fs, f0);
%
function [fh, nh] = CalculateHarmonics(signal, fs, f0)
    N = length(signal);
    Y = fft(signal, N);
    P2 = abs(Y / N);
    P1 = P2(1:floor(N / 2 + 1));

    f = fs * (0:(N / 2)) / N;

    [maxSpectrum, ~] = max(P1);
    maxSpectrum_dB = 20 * log10(maxSpectrum);

    threshold_dB = maxSpectrum_dB - 40;
    fmax_index = find(20 * log10(P1) >= threshold_dB, 1, 'last');
    fmax = f(fmax_index);

    nh = 0;
    fh = [];

    for n = 1:floor(fmax / f0)
        harmonicFreq = n * f0;
        [~, harmonicIndex] = min(abs(f - harmonicFreq));

        if f(harmonicIndex) <= fmax && P1(harmonicIndex) >= 0.1 * maxSpectrum
            nh = nh + 1;
            fh = [fh f(harmonicIndex)];
            fh = fh(end);
        end

    end

end

% AutocorrelationFundamentalFrequency calculates the fundamental frequency of a given signal using the autocorrelation method.
%
% Syntax:
%   f0_autocorr = AutocorrelationFundamentalFrequency(signal, fs)
%
% Input Arguments:
%   - signal: The input signal for which the fundamental frequency needs to be calculated.
%   - fs: The sampling frequency of the signal.
%
% Output Argument:
%   - f0_autocorr: The calculated fundamental frequency of the signal. If no peak is found in the autocorrelation, NaN is returned.
%
% Example:
%   signal = sin(2*pi*100*(0:1/fs:1));
%   fs = 1000;
%   f0 = AutocorrelationFundamentalFrequency(signal, fs);
%
% References:
%   [1] Rabiner, L. R., & Schafer, R. W. (1978). Digital processing of speech signals. Prentice-Hall, Inc.
function f0_autocorr = AutocorrelationFundamentalFrequency(signal, fs)
    [autocorr, lags] = xcorr(signal, 'coeff');
    % Find the first positive peak of the autocorrelation
    positiveLags = lags >= 0;
    positiveAutocorr = autocorr(positiveLags);
    positiveLags = lags(positiveLags);
    [pks, locs] = findpeaks(positiveAutocorr, 'MinPeakProminence', 0.3, 'MinPeakDistance', fs / 20);

    if ~isempty(locs)
        % Convert the lag to frequency
        peakLag = positiveLags(locs(1));
        f0_autocorr = fs / peakLag;
    else
        f0_autocorr = NaN;
    end

end


% GenerateFrequencyTables generates frequency tables for piano, violin, and flute.
% The function calculates the frequencies for each note in different octaves
% based on the base frequencies provided for each instrument.
%
% Output:
%   - tablePiano: Frequency table for piano notes in different octaves.
%   - tableViolin: Frequency table for violin notes in different octaves.
%   - tableFlute: Frequency table for flute notes in different octaves.
%
% Example usage:
%   [tablePiano, tableViolin, tableFlute] = GenerateFrequencyTables()
function [tablePiano, tableViolin, tableFlute] = GenerateFrequencyTables()
    notes = ["A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#"];

    octaves = ['0', '1', '2', '3', '4', '5', '6', '7'];

    baseFrequenciesPiano = [27.5, 29.135, 30.868, 16.35, 17.32, 18.35, 19.45, 20.60, 21.83, 23.12, 24.50, 25.96];
    baseFrequenciesViolin = [196, 207.65, 220, 233.08, 246.94, 130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00];
    baseFrequenciesFlute = [261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392, 415.3, 440, 466.16, 493.88];

    numOctaves = length(octaves);
    numNotes = length(notes);

    tablePiano = CalculateFrequencies(baseFrequenciesPiano, numOctaves, "piano");
    tableViolin = CalculateFrequencies(baseFrequenciesViolin, numOctaves, "violin");
    tableFlute = CalculateFrequencies(baseFrequenciesFlute, numOctaves, "flute");
end

function freqTable = CalculateFrequencies(baseFrequencies, numOctaves, type)
    freqTable = zeros(numOctaves, length(baseFrequencies));
    if type == "piano"
        for i = 1:numOctaves
            freqTable(i, :) = baseFrequencies * 2 ^ (i);
        end

    end

    if type == "violin"

        for i = 1:numOctaves
            freqTable(i, :) = baseFrequencies * 2 ^ (i - 3);
        end

    end

    if type == "flute"

        for i = 1:numOctaves
            freqTable(i, :) = baseFrequencies * 2 ^ (i - 4);
        end

    end

end

% CalculateSpectralCentroid calculates the spectral centroid of a given signal.
% The spectral centroid is a measure of the center of mass of the spectrum.
% 
% Inputs:
%   - signal: The input signal.
%   - fs: The sampling frequency of the signal.
%
% Output:
%   - spectralCentroid: The calculated spectral centroid.
%
% Example usage:
%   signal = [0.1, 0.2, 0.3, 0.4, 0.5];
%   fs = 44100;
%   centroid = CalculateSpectralCentroid(signal, fs);
function spectralCentroid = CalculateSpectralCentroid(signal, fs)
  spectrum = abs(fft(signal));
  normalizedSpectrum = spectrum / sum(spectrum);
  frequency = (0:length(spectrum) - 1) * (fs / length(spectrum));
  spectralCentroid = sum(frequency .* normalizedSpectrum);
end


% CalculateZeroCrossingRate calculates the zero crossing rate of a given signal.
% The zero crossing rate is defined as the number of times the signal changes sign divided by the length of the signal.
%
% Inputs:
%   - signal: The input signal for which the zero crossing rate needs to be calculated.
%
% Output:
%   - zeroCrossingRate: The calculated zero crossing rate of the signal.
%
% Example:
%   signal = [1, -2, 3, -4, 5];
%   zeroCrossingRate = CalculateZeroCrossingRate(signal);
function zeroCrossingRate = CalculateZeroCrossingRate(signal)
  zeroCrossingRate = sum(abs(diff(sign(signal)))) / length(signal);
end



% CalculateADSR calculates the Attack, Decay, Sustain, and Release parameters for a given signal.
% 
% Inputs:
%   - signal: The input signal
%   - fs: The sampling frequency of the signal
%   - maxSignal: The maximum value of the signal
% 
% Output:
%   - ADSRParams: A structure containing the calculated ADSR parameters
%
% Example usage:
%   signal = [0.1, 0.3, 0.5, 0.7, 0.9, 0.7, 0.5, 0.3, 0.1];
%   fs = 44100;
%   maxSignal = 1;
%   params = CalculateADSR(signal, fs, maxSignal);
% See also: DetectNoteTimes
function ADSRParams = CalculateADSR(signal, fs, maxSignal)
  % Calculate the start and end times of the note in the signal
  [startTemp, endTemp] = DetectNoteTimes(signal, fs, maxSignal);

  % Calculate the attack time
  peakAmplitude = max(signal);
  attackTime = find(signal >= peakAmplitude * 0.9, 1) / fs;

  % Calculate the decay time
  sustainLevel = peakAmplitude * 0.7;
  decayTime = find(signal <= sustainLevel, 1) / fs - attackTime;

  % Calculate the release time
  releaseStartIndex = floor(endTemp * fs);
  releaseTime = length(signal(releaseStartIndex:end)) / fs;

  % Store the calculated parameters in a structure
  ADSRParams.attack = attackTime;
  ADSRParams.decay = decayTime;
  ADSRParams.sustain = sustainLevel;
  ADSRParams.release = releaseTime;
end

function instrument = ClassifyInstrument(ADSRParams, spectralCentroid, zeroCrossingRate)
    % Define thresholds for classification based on empirical data
    spectralCentroidThresholds = struct('Piano', [500, 1500], 'Violin', [2000, 3000], 'Flute', [1500, 2500]);
    zeroCrossingRateThresholds = struct('Piano', [0.01, 0.02], 'Violin', [0.02, 0.05], 'Flute', [0.015, 0.025]);

    % Initialize instrument as 'Unknown'
    instrument = 'Unknown';

    % Check for Piano
    if (spectralCentroid >= spectralCentroidThresholds.Piano(1) && spectralCentroid <= spectralCentroidThresholds.Piano(2)) && ...
            (zeroCrossingRate >= zeroCrossingRateThresholds.Piano(1) && zeroCrossingRate <= zeroCrossingRateThresholds.Piano(2))
        instrument = 'Piano';
    end

    % Check for Violin
    if (spectralCentroid >= spectralCentroidThresholds.Violin(1) && spectralCentroid <= spectralCentroidThresholds.Violin(2)) && ...
            (zeroCrossingRate >= zeroCrossingRateThresholds.Violin(1) && zeroCrossingRate <= zeroCrossingRateThresholds.Violin(2))
        instrument = 'Violin';
    end

    % Check for Flute
    if (spectralCentroid >= spectralCentroidThresholds.Flute(1) && spectralCentroid <= spectralCentroidThresholds.Flute(2)) && ...
            (zeroCrossingRate >= zeroCrossingRateThresholds.Flute(1) && zeroCrossingRate <= zeroCrossingRateThresholds.Flute(2))
        instrument = 'Flute';
    end

end
