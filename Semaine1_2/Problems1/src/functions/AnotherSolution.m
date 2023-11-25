    
    % This algo is very simple, much simple than problemI and problem
    % I suggest that this algo is all you need to see, dont have to see the
    % other two funcitons it is too compicated(actually because i didnt 
    % really understand the question and i thought the question is
    % very difficult... thus i considered many boundary conditions
    % which is unnesscery)

    function AnotherSolution()
    % This function analyzes audio signals to identify intervals where the 
    % power exceeds a certain threshold. It visualizes these intervals on 
    % the signal and also plots the signal's power in dBm.

    % Add the path to the audio files
    addpath ../../../Audios

    % Read the audio files
    [s1, Fs1] = audioread("MarteauPiqueur01.mp3");
    [s2, Fs2] = audioread("Jardin01.mp3");
    [s3, Fs3] = audioread("Jardin02.mp3");
    [s4, Fs4] = audioread("Ville01.mp3");

    % Initialize parameters
    S = -48; % Sensitivity in dBV
    G = 40;  % Gain in dB
    P_SPL = 80; % Sound Pressure Level in dB SPL
    D_t = 1;   % Minimum duration threshold in seconds

    % Collect all signals and their sample rates
    signals = {s1, s2, s3, s4};
    Fs = {Fs1, Fs2, Fs3, Fs4};
    
    % Define window sizes for the analysis
    windowSizes = {39900, 31000, 40000, 37100};

    % Process each signal
    for i = 1:length(signals)
        % Calculate the minimum number of samples for the given D_t
        minSamples = D_t * Fs{i};

        % Initialize the array for storing cleaned intervals
        cleanedIntervals = []; 

        % Amplify the signal and convert SPL to dBm
        [P, interval] = SearchSignal(FunctionAmplifier(signals{i}, G), FunctionConvertSPLTodBM(P_SPL, S) + G, windowSizes{i});

        % Clean the intervals to ensure they meet the minimum duration
        for j = 1:size(interval, 1)
            duration = interval(j, 2) - interval(j, 1) + 1;
            if duration >= minSamples
                cleanedIntervals = [cleanedIntervals; interval(j, :)];
            end
        end

        % Compute the time vector for plotting
        t = (0:length(signals{i}) - 1) * (1 / Fs{i});  

        % Plot the original signal
        figure;
        subplot(2,1,1);
        plot(t, signals{i});
        hold on;

        % Mark the intervals that exceed the threshold
        halfWindowSize = floor(windowSizes{i} / 2);
        for k = 1:size(cleanedIntervals, 1)
            startIndex = max(1, (cleanedIntervals(k, 1) - halfWindowSize)); 
            endIndex = min(length(signals{i}), (cleanedIntervals(k, 2) + halfWindowSize)); 
            plot(t(startIndex:endIndex), signals{i}(startIndex:endIndex), 'r');
        end
        hold off;
        title('Signal with Over-Threshold Intervals Marked');
        xlabel('Time');
        ylabel('Signal Amplitude');

        % Plot the power in dBm
        subplot(2,1,2);
        numWindows = length(signals{i}) - windowSizes{i} + 1; 
        t = (0:numWindows - 1) * (1 / Fs{i});  
        p_dBm = CalculateWindowedPowerdBmSliding(signals{i}, windowSizes{i});
        plot(t, p_dBm);

        % Draw the threshold line
        yline(FunctionConvertSPLTodBM(P_SPL, S), 'r--');
    end
end

function p_dBm = CalculateWindowedPowerdBmSliding(signal, windowSize)
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

    % Convert the power from milliwatts to dBm
    p_dBm = 10 * log10(p_mW / 0.001);
end


function [P,interval] = SearchSignal(signal,seuilSpl, windowSize)
    % Search for intervals in the signal where the power exceeds a specified
    % threshold. Returns the power array and intervals exceeding the threshold.
    
    % Initialize variables
    interval = [];  % Array to store intervals that exceed the threshold
    penible = false;  % Flag to indicate whether currently in an interval
    lim = seuilSpl;   % Threshold in dBm

    % Calculate the power of the signal using sliding window
    P = CalculateWindowedPowerdBmSliding(signal, windowSize);
    begin_point = 0;

    % Iterate over the power array to find intervals exceeding the threshold
    for i = 1:length(P)
        if ~penible && P(i) >= lim
            penible = true;
            begin_point = i;  % Mark the start of an interval
        elseif penible && P(i) < lim
            penible = false;
            interval = [interval; [begin_point i]];  % Add the interval
        end
    end

    % Handle the case where the signal ends while still in an interval
    if penible
        interval = [interval; [begin_point length(P)]];
    end
end

function signal_amplified = FunctionAmplifier(signal, G)
    % Amplify the given signal by a specified gain.
    % The gain is provided in dB, and the amplification is applied linearly.
    
    % Convert gain from dB to linear scale and amplify the signal
    signal_amplified = signal * 10^(G / 20);
end


function P_dbm = FunctionConvertSPLTodBM(P_SPL, S)
    % Convert Sound Pressure Level (SPL) from dB to dBm.
    % This function takes SPL and sensitivity as inputs and returns the SPL in dBm.
    
    % Reference values for the conversion
    M_ref = 1;              % Reference magnitude
    P_reference = 20 * 10^(-6);  % Reference pressure in Pascals

    % Convert SPL to dBm
    P_dbm = 10 * log10((M_ref * 10^(S / 20) * P_reference * 10^(P_SPL / 20) * 10^2)^2 * 1000);
end
