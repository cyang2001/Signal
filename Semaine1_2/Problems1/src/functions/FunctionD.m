function [x, t] = FunctionD()
  F_e = 16000;
  f_0 = 1000;
  T_0 = 1/f_0;
  T_e = 1/F_e;
  phi = pi/3;
  p_dBm = 32;
  p_Watts = 1e-3 * 10^(p_dBm/10);
  B = sqrt(2*p_Watts);
  t = 0:T_e:5*T_0;
  x = B*sin(2*pi*f_0*t + phi);
  figure;
  plot(t, x);
  title("signal in the time domaine(x_5)");
  xlabel("second");
  ylabel("volt");
  frame = getframe(gcf);
  im = frame2im(frame);
  imwrite(im, '../../results/D.png')
end