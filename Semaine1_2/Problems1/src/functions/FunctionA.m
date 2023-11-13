function [x, t] = FunctionA()
    A = 2;
    f_0 = 1000;
    D = 2;
    F_e = 16000;
    T_e = 1/F_e;
    N = 5 * F_e / f_0;
    n = 0:N-1;
    t = n*T_e;
    x = A * sin(2*pi*f_0*n*T_e);
    figure;
    plot(n*T_e,x);
    title('Signal in the time domain - 5 periods');
    xlabel('seconds');
    ylabel('Volt');
end