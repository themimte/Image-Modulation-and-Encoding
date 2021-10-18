function [X,freq,spectrum]=centeredFFT(x,Fs) 
%this is a custom function that helps in plotting the two-sided spectrum 
%x is the signal that is to be transformed %Fs is the sampling rate 
 
N=length(x); 
 
if mod(N,2)==0     
         k=-N/2:N/2-1; % N even %this part of the code generates that frequency axis
else     k=-(N-1)/2:(N-1)/2; % N odd 
end
T=N/Fs; 
freq=k/T;  %the frequency axis
X=fft(x)/N; 
X=fftshift(X);
spectrum=abs(X).^2;
end