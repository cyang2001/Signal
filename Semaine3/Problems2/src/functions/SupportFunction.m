function SupportFunction()
    addpath ../../../Audios
    [s1, Fs1] = audioread('Pi_A_96K.wav');
    t = 0:1/Fs1:(length(s1)-1)/Fs1;
    s1 = mean(s1, 2);
    
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
