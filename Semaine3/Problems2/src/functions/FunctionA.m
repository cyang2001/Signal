function FunctionA()
  Fe = 8000;
  Te = 1/Fe;
  y = [];
  D = 8;%8ms
  t = 0:Te:D;
  for temp = t
    if((0<=temp) && (temp<=5))
      y = [y,1];
    else
      y = [y,0];
    end
  end
  figure;
  subplot(2,1,1);
  plot(t,y);
  title('y(n)')
  xlabel('ms');
  ylim([-0.5,1.5]);
  Y = fft(y);
  N = length(y);
  f = (-N/2:N/2-1)*(Fe/N);
  Y_shifted = fftshift(Y);
  subplot(2,1,2);
  plot(f,abs(Y_shifted));
  title('fft(f)')
  xlabel('Hz');
  ylabel('amplitude');
  frame = getframe(gcf);
  im = frame2im(frame);
  imwrite(im, '../../results/A.png');
end