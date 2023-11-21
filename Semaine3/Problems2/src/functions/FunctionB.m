function FunctionB()
  f1 = 1440;
  f2 = 2000;
  D = 1;
  Fe = 8000;
  Te = 1/Fe;
  t = 0:Te:D-Te;
  y1 = sin(2*pi*f1*t);
  y2 = sin(2*pi*f2*t);
  y1_2 = y1.*y2;
  Y = fft(y1_2);
  Y = fftshift(Y);
  N = length(Y);
  f = (-N/2:N/2-1)*Fe/N;
  figure;
  plot(f,abs(Y));
  xlabel('Frquence (Hz)');
  ylabel('Amplitude');
  
end