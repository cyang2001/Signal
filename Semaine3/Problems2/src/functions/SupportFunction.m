function SupportFunction()
    addpath ../../../Audios
    [s1, Fs1] = audioread('Pi_A_96K.wav');
    t = 0:1/Fs1:(length(s1)-1)/Fs1;
    plot(t, s1);
    xlabel('Time (s)');
    ylabel('Amplitude');
end 