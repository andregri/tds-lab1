%% Practica 1. Transfromada discreta de Fourier
% El objetivo de esta primera práctica es estudiar las propiedades de la 
% Transformada Discreta de Fourier y su inversa, así como claricar 
% conceptos como Resolución Espectral, Zero Padding, Enventanado, etc.
% Finalmente se aplicará el algoritmo de Goertzel para la detección de
% tonos en un sistema de alarmas.

%% Tarea 1
% $$ x(n) = a^{n} u(n) $$
%
% $$ x(z) = \frac{1}{1-az^{-1}} |a|<1 $$
%
% $$ z = e^{jw} $$
%
% $$ X(w) = \frac{1}{1-a e^{-jw}} $$

clear all;
close all;

a = input('Introduce a:');
N = input('Introduce N:');

n = 0:N-1;
x = a.^n;
wk = 2*pi*(0:N-1)/N;
X = 1./ (1-a*exp(-1i*wk));

xrec = ifft(X); 
xrecimag = imag(xrec);
xrecreal = real(xrec);

figure();
stem(n,xrecreal,'b');hold on;
stem(n,x,'ro');

hold off;
xlabel('n');
title(['a=',num2str(a),' N=',num2str(N)]);
legend('x(n)','x(n)_{rec}ifft')

% cuando el valor de a se acerca a 1, el error entre la señal recuperada 
% con la ifft y la señal original aumenta.

%% Tarea 2
clear all;
close all;

L = 5;
x = ones(1,L);
for N = [5 10 20 50 100 1000]
    n = 0:N-1;
    w = 2*n/N; % frecuencia normalizada en [0:1]
    X = abs(fft(x,N));
    figure;
    stem(w,X);
    xlabel('Normalised frequency');
    ylabel('Magnitude');
    title(['L=',num2str(L),' N=',num2str(N)]);
end

% El zero padding no añade información adicional, permite ver mejor el 
% perfil del espectro. Recupero siempre la misma señal.

%% Tarea 3
clear
close all
clc

Fm = 1000;
tt = [25 25 100 100];
N = [25 100 100 1000];

for i=1:4
    F1 = 90;
    F2 = 100;
    F3 = 240;
    F4 = 360;
    t = (Fm*tt(i))*0.001;
    x1 = cos(2*pi*(F1/Fm)*(0:t-1)) + cos(2*pi*(F2/Fm)*(0:t-1));
    x2 = x1 + cos(2*pi*(F3/Fm)*(0:t-1)) + cos(2*pi*(F4/Fm)*(0:t-1));
    y = fft(x2);
    % Generamos el eje X para interpretar frecuencias
    t1 = (0:t-1)*Fm/tt(i);
    figure;
    plot(t1,abs(y),'*-r');
    yy = fft(x2,N(i));
    hold on
    t2 = (0:N(i)-1)*Fm/N(i);
    plot(t2,abs(yy),'+-k')
    legend(['L=' num2str(t), ' N=' num2str(t)],['L=' num2str(t),' N=' num2str(N(i))])
    axis([0 Fm/2 0 max([abs(y) abs(yy)])]);
    grid
end
% t=25ms Fm=1000 L=25 N=25
% La resolución fisica es $$ R_f = \frac{Fm}{L} = 40Hz $$ pero la
% diferencia màs pequeña es de 10HZ. Por eso no podemos discernir las
% frecuencias F1=90Hz y F2=100Hz.

% t=25ms Fm=1000 L=25 N>25
% Aunque el orden de la DFT es mayor (la resolución computacional es mayor),
% la resolución fisica es siempre la misma. Entonces no podemos discernir 
% las frecuencias F1 y F2.

% t=100ms Fm=1000 L=100 N=100
% Hemos aumentado el numero de muestras de la señal L 4 veces y la
% resolución fisica es exactamente 10Hz. Ahora es posible ver F1 y F2

% t=100ms Fm=1000 L=100 N=1000
% si aumentamos N (zero padding), anadimos puntos en el espectro y aparecen 
% los puntos de la ventana rectangular superpuesta. Los picos en dos
% frecunecias diferentes de F1 y F2 porque los lobulos se solapan al espectro.

F = 10;
Fm = 1000;
tt = 250e-3;
L = Fm*tt;
N = L;
t = 0:L-1;
x = sin(2*pi*F/Fm*t);
subplot(211);
plot(t,x);
xlabel('t (ms)');
ylabel('x(t)');

X = fft(x,N);
w = Fm*(0:N-1)/N;
subplot(212);
stem(w,abs(X));
xlabel('f (Hz)');
ylabel('|X(w)|');

%% Tarea 4

clear
close all
N1 = 40; % Número de muestras de la señal
n1 = 0:N1-1;
F = 2*pi*[2000 2500 3000]/10000;
f1 = sum(cos(F'*n1));
N2 = 100; % Número de muestras de la señal
n2 = 0:N2-1;
f2 = sum(cos(F'*n2));
h_f1 = f1' .* blackman(N1);
h_f2 = f2' .* blackman(N2);
Nfft = 256;
ftf1 = fft(f1,Nfft); ftf2 = fft(f2,Nfft);
ftf3 = fft(h_f1,Nfft); ftf4 = fft(h_f2,Nfft);
w = [0:Nfft-1]*(2*pi/Nfft);
subplot(2,2,1);plot(w/pi,abs(ftf1),'k');xlabel('\omega / \pi');
ylabel('Magn. ');title('L=10, Rectangular)');
subplot(2,2,2);plot(w/pi,abs(ftf3),'k');xlabel('\omega / \pi');
ylabel('Magn. ');title('L=10, Hamming)');
subplot(2,2,3);plot(w/pi,abs(ftf2),'k');xlabel('\omega / \pi');
ylabel('Magn. ');title('L=20, Rectangular)');
subplot(2,2,4);plot(w/pi,abs(ftf4),'k');xlabel('\omega / \pi');
ylabel('Magn. ');title('L=20, Hamming)');

% Hay 3 frecuencias y la resta minima estre ellas es de 500Hz. Entonces
% hace falta una resolución fisica de 500Hz (necesitamos minimo L = 20).
% Sin embargo, usando la ventana de Hamming no podemos discernir las
% frecuencias cuando L = 20. Este efecto es debido al goteo espectral y a
% la anchura de la ventana. Además la amplitud de la señal enventanada con
% la ventana de Hamming es menor porque la frecuencias están distribuidas
% en un rango de frecuencias mayor.

% Gracias a las ventanas de Hamming, Hanning y Blackman el goteo espectral
% es menor porque la diferencia entre el primero y el segundo lobulo es
% mayor que la ventana rectangular. A pesar de un goteo espectral menor,
% para discernir frecuencias distintas, necesitamos una resolución fisica
% más grande por la anchura de la ventana.

%% Tarea 5
clear all;
close all;

N = 0:2048;
hold on;
plot(N,N.^2);
plot(N,N.*log2(N));
hold off;
% Hay una diferencia muy grande entre las funciones. La mejora es muy
% significativa.

x1 = [1 2 3];
x2 = [-1 1];
y = conv(x1,x2)
ydft1 = ifft(fft(x1,3) .* fft(x2,3), 3)
ydft2 = ifft(fft(x1,4) .* fft(x2,4), 4)
ydft3 = ifft(fft(x1,5) .* fft(x2,5), 5)
% Para calcular la convolución entre dos secuencias con la DFT, hemos usado
% esta formula: DFT = ifft(fft(x1,N) * fft(x2,N), N).
% Con N=3 no obtenemos el mismo resultado porque la longitud del resultado
% de la convolución es L = longitud_x1 + longitud_x2 - 1 = 4. Por eso
% tenemos aliasing en el primer termino.
% Con N=4 y con N=5 el resultado coincide con la convolución.

A = 0; B = 1; C = 0;
Fm = 500;
F1 = 150;
F2 = 175;
F3 = 200;
N = 20; % resolución computacional Rf = 25Hz
n = 0:N-1;
x = A*cos(2*pi*F1/Fm*n) + B*cos(2*pi*F2/Fm*n) + C*cos(2*pi*F3/Fm*n);

k1 = F1/Fm*N;
k2 = F2/Fm*N;
k3 = F3/Fm*N;

B1=[1];
A1=[1, -exp(i*2*pi*k1/N)];
y=filter(B1,A1,[x 0]);
y1=y(N+1)

B2=[1];
A2=[1, -exp(i*2*pi*k2/N)];
y=filter(B2,A2,[x 0]);
y2=y(N+1)

B3=[1];
A3=[1, -exp(i*2*pi*k3/N)];
y=filter(B3,A3,[x 0]);
y3=y(N+1)