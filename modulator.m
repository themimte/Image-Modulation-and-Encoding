function [moduled_pic_str,S] = modulator(pic_str,fc,fs,Ts,M)
num_bit=log2(M);
len=strlength(pic_str); % ttole ye string ba dastoore strlength gerefte mishe
if(mod(len,num_bit)~=0)
    for i=1:num_bit-mod(len,num_bit)
        pic_str=strcat(pic_str,"0");
    end    
end
len=strlength(pic_str); %toole string bad az inke be meghdare kafi sefr samte rast gharar dadim

num_of_t_periods=len/num_bit; 
moduled_pic_str=zeros(1,(Ts*fs+1)*num_of_t_periods); %khoroojie function
t=0:1/fs:Ts;

S=zeros(M,Ts*fs+1);
for m=1:M
    S(m,:)=sqrt(2/Ts)*cos(2*pi*fc.*t+ pi*(2*m-1)/M); %chon toole reshte kheli bolande va tedad pack ha ziade Sm ro tashkil midim
                                                    % ta har bar dige majboor nashim Sm ro vase har packhesab konim
end
%collect=strings([num_of_t_periods,1]); % TEST
%collect=zeros(num_of_t_periods,1); % TEST
%k=1;   % TEST
temp=strings([1,1]); %vaghti mikhaym ye araye m*n tarif konim ke har derayash string bashe az dastoor strings([m,n]) estefade mikonim

down_interval=zeros(1,1); %hade payeene bazehaye zamani motanaseb ba modulation har packet
up_interval=zeros(1,1);   %hade balaye bazehaye zamani motanaseb ba modulation har packet 
for i=1:num_bit:len
    down_interval=[down_interval,1+((Ts*fs+1)*(i-1))/num_bit]; %concatenate
    up_interval=[up_interval,((Ts*fs+1)*(i-1+num_bit))/num_bit];    %concatenate
end
down_interval=down_interval(1,2:end);   %avalin hade payeen 0 hast ke ma nemikhaym va baghie catenate shode haro mikhaym pas az 2vomi barmidarim
up_interval=up_interval(1,2:end);   %avalin hade bala 0 hast ke ma nemikhaym va baghie catenate shode haro mikhaym az 2vomi barmidarim

for i=1:num_of_t_periods
    temp =pic_str{1,1}(1+(i-1)*num_bit:i*num_bit); %harbar pack be toole num_bit ro barmidarim va moduleh mikonim
    %collect{k,1}=temp;k=k+1; % TEST
    binar=str2num(temp); m_prime=bin2dec(string(binar)); %moadele decimal packi ke bardashtim
   %collect(k,1)=m_prime;k=k+1; % TEST
    moduled_pic_str(1,down_interval(i):up_interval(i))=S(m_prime+1,:); %bejaye oon pack meghdare moduleh shode ro too moduled_pic_str mizarim
end

end