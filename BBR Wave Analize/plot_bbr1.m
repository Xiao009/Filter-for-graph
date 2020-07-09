#clean up history and variable 
clc
clear all

#Open upn file as array
[fname, fpath, fltidx] = uigetfile("*.log");
Fs = 214.65;

pkg load signal
[time force deflection] = textread(strcat(fpath,fname),"%d %d %d");
if (time(2)==time(3))
  buffer = time;
  time = force;
  force = buffer;
endif
force = force(2:end);
length(force)
count=0;
first=0;
x1=0;
x2=0;
for i=1:7
if (first==1)
  if (force(i)==force(i+1))
    count++;
   else
    count++;
    x2=i;
    break;
  endif
else
   if (force(i)!=force(i+1))
  count=0;
  first=1;
  x1=i;
  endif
endif
endfor
T = (time(x2)-time(x1));
Fs = double(200000.0)/double(T) ;
dt=0:T:5e-3-T;
load = downsample(force,count);
x=load
subplot(2,1,1);
plot(double(x)/double(10000.0));
title('initial signal');

#fft, frequency domain transfer
k1 = load([200:1100]);
L=length(k1)
y=fft(k1)
P2 = abs(y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

yfft=fftshift(fft(y));
subplot(2,1,2);
plot(f,P1);
axis([0 50 0 100]);
title('frequence');
set(gca,'XTick',0:2:100)

function my_output=ideallp(wc,N) 

%Ideal Lowpass filter computation 

% hd=ideal impulse response between 0 to M-1 

% wc=cutoff frequency in radians 

% M=length of the ideal filter 


alpha=(N-1)/2; 

n=0:1:(N-1); 

m=n-alpha+eps; 

my_output=sin(wc*m)./(pi*m); 

end 


#filter used
#Fs is around 4.4866 sps
wp=0.5;       # passband  fs=0.5              
ws=0.64;       # stopband   fp=1             
deltaw=ws-wp;    #transition zone width                 
N=ceil(8.9*pi/deltaw);  #1.8  max 8.9   ceil(8.9*pi/deltaw)     
wc=(ws+wp)/2          #bandpass range
hd=ideallp(wc,N);                 
wd1=blackman(N);   #23db, so pick blackman(25db)
wd2=hamming(N);
h1=hd.*wd1';
h2=hd.*wd2';
#final=fftfilt(h1,double(x)/double(10000.0),200);
final=conv(h1,double(k1)/double(10000.0));
[H1,w]=freqz(h1,1);
[H2,w]=freqz(k1,1);    
[H4,w]=freqz(h2,1);

result=conv(H1,H2);
result1=conv(H1,H4);
[H3,w]=freqz(result,1);
[H5,w]=freqz(result1,1);  


figure;
subplot(2,1,1);
plot(final,'r',double(x)/double(10000.0),'b');
legend('processed signal','original signal');
title('Signal Magnitude');
subplot(2,1,2);
plot(w/pi,20*log10(abs(H3))-220,'r',w/pi,20*log10(abs(H2)),'b');
legend('blackman','test object');
title('frequence response');


# kalman method
k=double(load)/double(10000.0);
# calculate varriance
mean=sum(k)/1149
Num=1149
sum_diff=sum((k-mean).^2)
variance=sum_diff/Num;    #2
Value_ideal=99.9;
size=[Num,1];
value_Start=0;
Q=1e-3;  #1e-3  trust the ideal value 
R=0.36; #0.36 not trust the sensor value
M_kalman(74)=0;    # set up the start point to start calculation for Calrman
P_kalman(74)=variance;
array_ideal=Value_ideal*ones(size);
#array_ideal=array_ideal([74:1149]);

#Kalman: Best evaluation
#Pre: prediction
for G=75:Num
  M_pre(G)=M_kalman(G-1);
  P_pre(G)=P_kalman(G-1)+Q;
  K(G)=P_pre(G)/(P_pre(G)+R);
  M_kalman(G)=M_pre(G)+K(G)*(k(G)-M_kalman(G-1));
  P_kalman(G)=P_pre(G)-K(G)*P_pre(G);
end
plot(Value_ideal*ones(size),'g',k,'b',M_kalman,'r');



#ynfft1=fftshift(fft(final)) 
#ynfft1=bandstop(ynfft1,[600 615])
#figure;
#subplot(2,1,1);
#plot(abs(ynfft1))
#title('frequence filted');
#orig=abs(y)
#subplot(2,1,2)
#plot(orig)

