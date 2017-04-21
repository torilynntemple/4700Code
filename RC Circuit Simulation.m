%fdghfd

clear 
clc


%% QUESTION 1

R1 = 20;
C1 = 10*10^-6; 
wo = 1/(R1*C1);
fc = 1/(2*pi*R1*C1); 

figure (1)
f = linspace(0,10000,1000);
transfer_fn = 1./(abs((1j.*2.*pi.*f.*R1.*C1)+1)); 
logtransfer = 20.*log10(transfer_fn);
figure (1)
semilogx(f,logtransfer)
title('Transfer Function Outlining the Cutoff Frequency of the Low Pass Filter');
xlabel('Frequency (Hz)');
ylabel('Transfer Function (dB)');

V1 = 1; 
t = linspace(0,0.1, 1000);
Vo = V1.*(1-exp(-1.*t./(R1.*C1)));

%for question blahblah
figure(2)
plot(t,Vo)
title('Output Voltage With Respect to Time of a Low Pass Filter');
xlabel('Time (s)');
ylabel('Voltage (V)');
hold on

%for question 1g)
%when testing the code below for varying time steps, only change N. 
Vi=1;
N = 100; %numbe of time step
timestep = linspace(0,0.002,N);
time = 0.002/N; 
Voltage1 = zeros(N,1); 

    for i = 1:(N-1)
        Voltage1(i+1) = Voltage1(i) + time*(Vi-Voltage1(i))/(R1*C1);
    end
    
figure (3)
plot(timestep,Voltage1)
hold on

%for question 1h)
N = 1000; %numbe of time step
timestep = linspace(0,0.002,N);
f = 6000; 
time = 0.002/N; 
Voltage = zeros(N,1); 

    for i = 1:(N-1)
        Voltage(i+1) = Voltage(i) + time*(sin(2*pi*f*time*i)-Voltage(i))/(R1*C1);
    end
    
figure (4)
plot(timestep,Voltage)
title('Output Voltage With Respect to Time of a Low Pass Filter');
xlabel('Time (s)');
ylabel('Voltage (V)');

% 
% %for question 1i)
figure (5)
y = abs(fftshift(fft(Voltage)));
freq_axis = (-N/2:N/2-1)/(time*N);
%freq = linspace(-7500,7500,N)
plot (freq_axis, y)
title('Output Voltage With Respect to Time of a Low Pass Filter');
xlabel('Frequency (Hz)');
ylabel('DFT');

%% QUESTION 3
V1 = 1; 
R1 = 20;
C1 = 10*10^-6; 
wo = 1/(R1*C1);
fc = 1/(2*pi*R1*C1); 
Zc1 = 1/(j*2*pi*fc*C1); 
iterations = 1000; 
Imax = Vi/(5*R1); 
timestep = linspace(0,0.002, iterations);
time = 0.002/iterations; 
Imax = Vi/(5*R1); 
Vmax = Imax*R1; 
Vrand = randn(1,length(timestep));
Vrand = Vrand.*Vmax./max(Vrand);  

Vo = (V1+ Vrand).*(1-exp((-1.*timestep)./(R1.*C1)));

figure (6)
plot(timestep, Vo) 
title('Voltage Output With the Addition of Current Noise ');
xlabel('Time (s)');
ylabel('Voltage (V)');
 
%for question 3b)

RMS_Vo =sqrt( mean(Vo.^2))


%for question 3c)
figure (7)
y = abs(fftshift(fft(Vo)));
freq_axis = (-iterations/2:iterations/2-1)/(time*iterations);
plot (freq_axis, y)
title('Output Voltage With Respect to Time of a Low Pass Filter with Noise');
xlabel('Frequency (Hz)');
ylabel('DFT');


%for question 3d
V1 = 1; 
R1 = 10;
C1 = 10*10^-6; 
iterations = 1000; 
Imax = Vi/(5*R1); 
timestep = linspace(0,0.002, iterations);
time = 0.002/iterations; 
Imax = Vi/(5*R1); 
Vmax = Imax*R1; 
Vrand = randn(1,length(timestep));
Vrand = Vrand.*Vmax./max(Vrand);  

Vo = (V1+ Vrand).*(1-exp((-1.*timestep)./(R1.*C1)));
RMS_Vo =sqrt( mean(Vo.^2))

figure (9)
y = abs(fftshift(fft(Vo)));
freq_axis = (-iterations/2:iterations/2-1)/(time*iterations);
%freq = linspace(-7500,7500,N)
plot (freq_axis, y)
title('Frequency Spectrum with R = 10 ohms and C = 10uF');
xlabel('Frequency (Hz)');
ylabel('DFT');


R = 50;
figure (10)
plot (freq_axis, y)
title('Frequency Spectrum with R = 50 ohms and C = 10uF');
xlabel('Frequency (Hz)');
ylabel('DFT ');

R = 20; 
C = 20*10^-6;
figure (11)
plot (freq_axis, y)
title('Frequency Spectrum with R = 20 ohms and C = 20uF');
xlabel('Frequency (Hz)');
ylabel('DFT');
