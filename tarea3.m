function tarea3(Fm,tt,N)
    F1 = 90;
    F2 = 100;
    F3 = 240;
    F4 = 360;
    t = (Fm*tt)*0.001;
    x1 = cos(2*pi*(F1/Fm)*(0:t-1)) + cos(2*pi*(F2/Fm)*(0:t-1));
    x2 = x1 + cos(2*pi*(F3/Fm)*(0:t-1)) + cos(2*pi*(F4/Fm)*(0:t-1));
    y = fft(x2);
    % Generamos el eje X para interpretar frecuencias
    t1 = (0:t-1)*Fm/tt;
    plot(t1,abs(y),'*-r');
    yy = fft(x2,N);
    hold on
    t2 = (0:N-1)*Fm/N;
    plot(t2,abs(yy),'+-k')
    legend(['L=' num2str(t), ' N=' num2str(t)],['L=' num2str(t),' N=' num2str(N)])
    axis([0 Fm/2 0 max([abs(y) abs(yy)])]);
    grid
end