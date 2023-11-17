function FunctionC()
  [x,t] = FunctionA;
  x_3 = abs(x);
  x_4 = [];
  T_0 = 1/1000;
  T_e = 1/16000;
  t_ = 0:T_e:0.005;
  for temp=t_
    if(mod(temp,T_0) >= 0 && mod(temp,T_0) < (T_0/2))
      x_4=[x_4,2];
    else
      x_4=[x_4,0.5];
    end
  end
  figure;
  subplot(2,1,1);
  plot(t,x_3);
  title("Signal in the time domain(x_3)");
  xlabel("second");
  ylabel("volt");
  subplot(2,1,2);
  plot(t_,x_4);
  title("Signal in the time domain(x_4)");
  xlabel("second");
  ylabel("volt");
end