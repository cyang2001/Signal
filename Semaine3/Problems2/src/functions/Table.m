function Table()

    octaves = ['0', '1', '2', '3', '4', '5', '6', '7'];
    notesPiano = ["A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#"];
    notesViolin = ["G", "G#", "A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#"];
    notesFlute = ["C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#", "A", "Bb", "B"];
    [tablePiano, tableViolin, tableFlute] = generateFrequencyTables();


    for i = 1 : length(octaves)
      for j = 1 : length(notesViolin)
        disp(["octave: ", num2str(octaves(i)) ,"note: ", notesViolin(j) , num2str(tableViolin(i, j))])

      end

    end
  end
  function [tablePiano, tableViolin, tableFlute] = generateFrequencyTables()
    notes = ["A", "Bb", "B", "C", "C#", "D", "Eb", "E", "F", "F#", "G", "G#"];

    octaves = ['0', '1', '2', '3', '4', '5', '6', '7'];

    baseFrequenciesPiano = [27.5, 29.135, 30.868, 32.703, 34.648, 36.708, 38.891, 41.203, 43.654, 46.249, 48.999, 51.913];
    baseFrequenciesViolin = [196, 207.65, 220, 233.08, 246.94, 130.81, 138.59, 146.83, 155.56, 164.81, 174.61, 185.00];

    baseFrequenciesFlute = [261.63, 277.18, 293.66, 311.13, 329.63, 349.23, 369.99, 392, 415.3, 440, 466.16, 493.88];

    numOctaves = length(octaves);
    numNotes = length(notes);

    tablePiano = calculateFrequencies(baseFrequenciesPiano, numOctaves, "piano");
    tableViolin = calculateFrequencies(baseFrequenciesViolin, numOctaves, "violin");
    tableFlute = calculateFrequencies(baseFrequenciesFlute, numOctaves, "flute");
end

function freqTable = calculateFrequencies(baseFrequencies, numOctaves, type)
    freqTable = zeros(numOctaves, length(baseFrequencies));
    if type == "piano"
      for i = 1:numOctaves
        freqTable(i, :) = baseFrequencies * 2^(i);
      end
    end
    if type == "violin"
      for i = 1:numOctaves
        freqTable(i, :) = baseFrequencies * 2^(i - 3);
      end
    end
    if type == "flute"
      for i = 1:numOctaves
        freqTable(i, :) = baseFrequencies * 2^(i - 4);
      end
    end
end

