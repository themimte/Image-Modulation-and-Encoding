function [demodulated_picture] = demodulation(modulated_picture,fc,fs,Ts,M)
%input is modulated signal
%output is demodualted signal i.e. seperated signal
%using matched filter
rroww=length(modulated_picture)/(Ts*fs+1);
demodulated_picture=zeros(rroww,Ts*fs+1);
for i=1:round(rroww)
          demodulated_picture(i,:)=modulated_picture( i*(Ts*fs+1)-Ts*fs : i*(Ts*fs+1) );
end
end

