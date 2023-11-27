function FunctionA()
  Fe = 8000;
  Te = 1/Fe;
  y = [];
  D = 0.006;
  n = 1 : 480;
  t = Te*n;
  for temp = t
    if((0<=temp) && (temp<=0.005))
      y = [y 1];
    else
      y = [y 0];
    end
  end
  figure;
  subplot(2,1,1);
  plot(t,y);
  title('y(n)')
  xlabel('s');
  Y = fft(y);
  N = length(y);
  f = (0:N-1)*(Fe/N);
  %Y_shifted = fftshift(Y);
  subplot(2,1,2);
  plot(f,abs(Y));
  xlim([0 4000]);
  title('fft(f)')
  xlabel('Hz');
  ylabel('amplitude');
  frame = getframe(gcf);
  im = frame2im(frame);
  imwrite(im, '../../results/A.png');
end