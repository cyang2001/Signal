function [s1,Ts1,s2,Ts2,s3,Ts3,s4,Ts4] = FunctionF()
  addpath ../../../Audios
  [s1, Fs1] = audioread("MarteauPiqueur01.mp3");
  [s2, Fs2] = audioread("Jardin01.mp3");
  [s3, Fs3] = audioread("Jardin02.mp3");
  [s4, Fs4] = audioread("Ville01.mp3");
  n1 = length(s1);
  Ts1 = 1/Fs1;
  t1 = 0:Ts1:(n1-1)*Ts1;
  n2 = length(s2);
  Ts2 = 1/Fs2;
  t2 = 0:Ts2:(n2-1)*Ts2;
  Ts3 = 1/Fs3;
  n3 = length(s3);
  t3 = 0:Ts3:(n3-1)*Ts3;
  n4 = length(s4);
  Ts4 = 1/Fs4;
  t4 = 0:Ts4:(n4-1)*Ts4;
  Fenetre(s1, 10000);
  power_s1 = mean(s1.^2);
  power_s2 = mean(s2.^2);
  power_s3 = mean(s3.^2);
  power_s4 = mean(s4.^2);
  %fprintf('power of s1 is: %f w\n', power_s1);
  %fprintf('power of s2 is: %f w\n', power_s2);
  %fprintf('power of s3 is: %f w\n', power_s3);
  %fprintf('power of s4 is: %f w\n', power_s4);
  figure;
  subplot(4,1,1);
  s1_dBm = FunctionCalculerPowerMeandBM(s1)
  plot(t1,s1);
  title('signal of MarteauPiqueur01')
  xlabel('second');
  ylabel('v');
  s2_dBm = FunctionCalculerPowerMeandBM(s2)
  subplot(4,1,2);
  plot(t2,s2);
  title('signal of Jardin01');
  xlabel('second');
  ylabel('v');
  s3_dBm = FunctionCalculerPowerMeandBM(s3);
  subplot(4,1,3);
  plot(t3,s3);
  title('signal of Jardin02');
  xlabel('second');
  ylabel('v');
  s4_dBm = FunctionCalculerPowerMeandBM(s4);
  subplot(4,1,4);
  plot(t4,s4);
  title('signal of Vill01');
  xlabel('second');
  ylabel('v');
  frame = getframe(gcf);
  im = frame2im(frame);
  %imwrite(im, '../../results/F.png');
figure;
windowSize = 1024;
p_dBm = calculateWindowedPowerdBm(s1, windowSize);

t = (0:(length(p_dBm) - 1)) * (windowSize / Fs1);

plot(t, p_dBm);
xlabel('Time (s)');
ylabel('Power (dBm)');
title('Signal Power in dBm');
hold on;
yline(8, 'r--');
hold off;
grid on;

end

function power = Fenetre(signal, lengthF)
  N = length(signal);
  power = []; 
  for i = 1:lengthF:N-lengthF+1
    window = signal(i:i+lengthF-1);  
    windowPower = mean(window.^2);  
    power = [power, windowPower];   
  end
  for j = 1:length(power)
    %fprintf('power of window%d is %f\n',j,power(j));
  end
  %fprintf('mean of power is %f\n',mean(power));
end

function power_mean_dBm = FunctionCalculerPowerMeandBM(signal)
  power_mean_dBm = 10 * log10(FunctionCalculerPowerMeanmW(signal) / 0.001);
end

function power_mean_mW = FunctionCalculerPowerMeanmW(signal)
  power_mean_mW = mean(signal .^ 2) / 1000;
end
function p_dBm = calculateWindowedPowerdBm(signal, windowSize)
    numWindows = floor(length(signal) / windowSize);
    p_mW = zeros(1, numWindows);
    
    for i = 1:numWindows
        windowStart = (i - 1) * windowSize + 1;
        windowEnd = i * windowSize;
        window = signal(windowStart:windowEnd);
        p_mW(i) = mean(window.^2);
    end
    p_dBm = 10 * log10(p_mW / 0.001);
end

