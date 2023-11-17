function [y, t2] = FunctionA()
    A = 2;
    f_0 = 1000;
    D = 2;
    F_e = 16000;
    T_e = 1/F_e;
    t_0 = 1/f_0;
    t = 0:T_e:D;
    p1 = 0;
    p2 = 5;
    [~,i1] = min(abs(t-p1*t_0));
    [~,i2] = min(abs(t-p2*t_0));
    x = A * sin(2*pi*f_0*t);
    t2 = t(i1:i2);
    y = x(i1:i2);
    figure;
    plot(t2,y);
    title('Signal in the time domain - 5 periods');
    xlabel('seconds');
    ylabel('Volt');
    ylim([-2.25,2.25]);
    grid on; 
end