function FunctionA()
    A = 2;
    f_0 = 1000;
    D = 2;
    F_e = 16000;
    T_e = 1/F_e;
    n = 0:T_e:5*1/f_0;
    x = A * sin(2*pi*f_0*n);
    figure;
    plot(n,x);
    title('Signal in the time domain - 5 periods');
    xlabel('seconds');
    ylabel('Volt');