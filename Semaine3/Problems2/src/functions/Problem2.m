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
    audioName = {'Pi_A_96K.wav', 'Pi_C_96K.wav', 'Vi_A3_96K.wav', 'Vi_C3_96K.wav', 'Vi_G4_96K.wav', 'Fl_A4_96K.wav', 'Fl_B3_96K.wav'};
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
        
        % Detect note begin and end times
        [startTemp, endTemp] = DetectNoteTimes(p_W, Fs{i}, max(p_W));

        disp(['The begining of this signal is :', num2str(startTemp), ' ending is :', num2str(endTemp)]);

        % Calculate average power  P
        avgPower = CalculateAveragePower(audioSignal{i}, Fs{i}, startTemp, endTemp);

        disp(['The average power of this note is :', num2str(avgPower), 'W']);

        % Calculate fundamental frequency
        f0 = AutocorrelationFundamentalFrequency(audioSignal{i}, Fs{i});

        disp(['The fundamental frequency of this note is :', num2str(f0), 'Hz']);
        % Only for test
        if i == 1 || i == 2
            tolerance = 2;
            [octave, note] = ind2sub(size(tablePiano), find(abs(tablePiano - f0) <= tolerance));
            disp(['The note of this frequency is :', notesPiano(note)]);
            disp(['The octave of this frequency is :', octaves(octave)]);
        end

        if i == 3 || i == 4 || i == 5
            tolerance = 4;
            [octave, note] = ind2sub(size(tableViolin), find(abs(tableViolin - f0) <= tolerance));
            disp(['The note of this frequency is :', notesViolin(note)]);
            disp(['The octave of this frequency is :', octaves(octave)]);
        end

        if i == 6 || i == 7
            tolerance = 2;
            [octave, note] = ind2sub(size(tableFlute), find(abs(tableFlute - f0) <= tolerance));
            disp(['The note of this frequency is :', notesFlute(note)]);
            disp(['The octave of this frequency is :', octaves(octave)]);

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
        ylabel('Power/Frequency (dB/Hz)');
        title('Power Spectral Density of the Signal ');
        frame = getframe(gcf);
        im = frame2im(frame);

        imwrite(im, ['../../results/' , audioName{i} , '.png']);
        %instrument = ClassifyInstrument(CalculateADSR(p_W, Fs{i}, max(p_W)), CalculateSpectralCentroid(p_W, Fs{i}), CalculateZeroCrossingRate(p_W));
        %disp([instrument]);
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

function avgPower = CalculateAveragePower(audioSignal, fs, startTime, endTime)
    startIndex = floor(startTime * fs);
    endIndex = ceil(endTime * fs);
    avgPower = mean(audioSignal(startIndex:endIndex) .^ 2) / 0.001;

end

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

function f0_autocorr = AutocorrelationFundamentalFrequency(signal, fs)
    % Calculate the fundamental frequency of the given signal using the
    % autocorrelation method
    % Calculate the autocorrelation of the signal
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

function spectralCentroid = CalculateSpectralCentroid(signal, fs)
    spectrum = abs(fft(signal));
    normalizedSpectrum = spectrum / sum(spectrum);
    frequency = (0:length(spectrum) - 1) * (fs / length(spectrum));
    spectralCentroid = sum(frequency .* normalizedSpectrum);
end

function zeroCrossingRate = CalculateZeroCrossingRate(signal)
    zeroCrossingRate = sum(abs(diff(sign(signal)))) / length(signal);
end

function ADSRParams = CalculateADSR(signal, fs, maxSignal)

    [startTemp, endTemp] = DetectNoteTimes(signal, fs, maxSignal);

    % Attack time: Time taken for the amplitude to reach its peak after onset
    peakAmplitude = max(signal);
    attackTime = find(signal >= peakAmplitude * 0.9, 1) / fs;

    % Decay time: Time taken for the amplitude to drop to sustain level
    sustainLevel = peakAmplitude * 0.7; % Assume sustain level is 70 % of peak
    decayTime = find(signal <= sustainLevel, 1) / fs - attackTime;

    % Release time: Time taken for the sound to fade away after the end of the note
    % For simplicity, we can assume that the release phase starts at the end of the note
    releaseStartIndex = floor(endTemp * fs);
    releaseTime = length(signal(releaseStartIndex:end)) / fs;

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
