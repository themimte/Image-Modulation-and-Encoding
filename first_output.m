%phase 1 source encoder,decoder
%% algorithm
 clear all;close all;clc;
  a=imread('36.gif');
  a=imresize(a,0.125);
  a=double(a);
%a=[1,1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4;5,5,6,6,7,7,8,8,9,9,10,10,11,12,13,14,15,16];
[C,ia,ic]=unique(a);
a_counts = accumarray(ic,1);
value_counts = [C, a_counts];
sorted_value_counts = sortrows(value_counts, -2);
final_sorted_value_counts=sorted_value_counts; %used in the final level of encoding
grayscales=sorted_value_counts(:,1);
sorted_value_counts=sorted_value_counts(:,2);
temp=sorted_value_counts;
[row , column] = size(sorted_value_counts);
min_freq=zeros(row-1,3);%[sum of 2 mins;big min;little min]

for i=1:row-2
    min_freq(i+1,1)=sorted_value_counts(row-i) + sorted_value_counts (row-i+1);
    min_freq(i+1,2)=max(sorted_value_counts(row-i),sorted_value_counts (row-i+1));
    min_freq(i+1,3)=min(sorted_value_counts(row-i),sorted_value_counts (row-i+1));
    sorted_value_counts=sorted_value_counts(1:row-i-1,:);
    sorted_value_counts=cat(1,sorted_value_counts,min_freq(i+1));
    sorted_value_counts=sort(sorted_value_counts,'descend');  
end
min_freq=sort(min_freq,'descend');
min_freq=cat(1,[sum(sorted_value_counts),max(sorted_value_counts),min(sorted_value_counts)] , min_freq(1:row-2,:) );

coding=zeros(row-1,3);coding=coding+3;coding(1,1)=0;coding(1,2)=1;coding(1,3)=0;
coding=string(coding); coding{1,1}='-';
for i=1:row-2
    for j=i+1:row-1
        
        if coding{j,1} == '3' & min_freq(i,2)== min_freq(j,1)  
          
                coding{j,1} = coding{i,2} ;
                coding{j,2} = [coding{j,1},'1'];
                coding{j,3} = [coding{j,1},'0'];
                break;
        end
        
    end
      for j=i+1:row-1
        if  coding{j,1}=='3' & min_freq(i,3)==min_freq(j,1)
                coding{j,1}=coding{i,3};
                coding{j,2}=[coding{j,1},'1'];
                coding{j,3}=[coding{j,1},'0'];
               break;       
        end 
       end
end

alter=coding;
alter=string(alter);
k=1;flag1=1  ;  flagig1=1 ; flag2=1 ; flagig2=1;
for i=1:row-2
    flag1=1;
    flag2=1;
    flagig1=1;
    flagig2=1;
    for j=i+1:row-1
        
        if (flagig1==1)  & ( size(alter{i,2}) == size( alter{j,1} ) )
           if alter{i,2}==alter{j,1}
                flag1=0;
                flagig1=0;
                alter(j,1)='3';
                alter(i,2)='3';
           end
        end
        if (flagig2==1) & ( size(alter{i,3}) == size( alter{j,1} ) )
            if alter{i,3}==alter{j,1}
                flagig2=0;
                flag2=0;
                alter(j,1)='3';
                alter(i,3)='3';
            end
        end
    end
end
k=1;
code=zeros(row,1);code=code+3;code=string(code);
for i=1:row-1
    if alter{i,2}~='3'
        code{k,1}=alter{i,2};
        k=k+1;
    end
        if alter{i,3}~='3'
         code{k,1}=alter{i,3};
         k=k+1;
        end    
end
str_value_counts=string(final_sorted_value_counts);
str_value_counts=cat(2,str_value_counts,code);
name=["grayscale","freq","code"];
str_value_counts=cat(1,name,str_value_counts);

%%encoder
pic_code=zeros(size(a)); pic_code=string(pic_code);
[row_num,~]=size(str_value_counts); row_num=row_num-1;
[r,cc]=size(a);
for i=1:r
    for j=1:cc
        for k=2:row_num+1
            if a(i,j)==str2num(str_value_counts{k,1})
                pic_code{i,j}=str_value_counts{k,3};
            end
        end
    end
end

%% Modulator
%pic_code=["1010110","111001","10101010";"000101011","00110101","11001010010100"];
pic_str=join(join(pic_code)'); pic_str=erase(pic_str," "); %all  elements of pic_code comes together in a row vector with name 'pic_str'
len_pic_str=strlength(pic_str);
fs=100000;Ts=.01;fc=10000;M=4;
[modulated_picture,Sm]=modulator(pic_str,fc,fs,Ts,M);
t=0:1/fs:Ts;
%plot(t,modulated_picture(1,1:length(t)));
%% demodulator 
demodulated_picture=demodulation(modulated_picture,fc,fs,Ts,M);
%% detector
rcvd_symbols=detector(demodulated_picture,Sm,fc,fs,Ts,M);

rcvd_finaly_string=join(join(rcvd_symbols)'); rcvd_finaly_string=erase(rcvd_finaly_string," ");
rcvd_finaly_string=rcvd_finaly_string{1,1}(1:len_pic_str);

namoosan_kar_kon=strcmp(rcvd_finaly_string,pic_str);


%% decoder
size_a=size(a);
final_pic_array = decoder(rcvd_finaly_string,code,grayscales,size_a);

in_tan_bemire_kar_kon=(max(max(final_pic_array-a))==0);
 
  final_pic_array=uint8(final_pic_array);     %felan niazi besh nadarim
  final_pic_array=imresize(final_pic_array,8);    %felan niazi besh 
  figure
  imshow(final_pic_array);   %felan niazi besh nadarim


%% phase 2
%% part B
[X,freq,spectrum]=FFT_function(modulated_picture,fs);
%% Part C
energy=sum(spectrum);
figure
plot(freq,abs(X));
title('frequency spectrum of 1.gif');
FC=find(spectrum==max(spectrum));
FC_pos=FC(1,2);
FC_neg=FC(1,1);
B_W_energy=spectrum(FC_pos);
for i=1:length(freq)-FC_pos
    if B_W_energy>.99*energy
        break;
    else
        B_W_energy = B_W_energy + 2*( spectrum(FC_pos+i)+spectrum(FC_pos-i) );
    end
end
B_W=i*fs/length(freq);
FC_channel=freq(FC_pos);
fprintf('for 1.gif fc of band-pass filter is %d Hz and BW is %d\n',FC_channel,B_W);
frequency=cat(2, freq(FC_neg-i:FC_neg+i ),zeros(1,2*i-1),freq(FC_pos-i:FC_pos+i) );
spec=cat(2, spectrum(FC_neg-i:FC_neg+i ),zeros(1,2*i-1),spectrum(FC_pos-i:FC_pos+i) );
figure
plot( frequency ,spec);
title('Spectrum density for filtered 1.gif');
ylabel('S(x)');xlabel('f');
%% Part D

%modulated_noisy_picture=normrnd(0,10,size(modulated_picture)) + modulated_picture;
modulated_noisy_picture=modulated_picture;%zeros(length(modulated_picture),1);
LLL=length(modulated_picture);
varia=linspace(1,1000,40);
error_probability=zeros(1,40);
last_len=length(rcvd_finaly_string);
varian=40;error_prob=0;
%for cntr=1:40
%noise=normrnd(0,varia(cntr),[1,LLL]);
noise=normrnd(0,varian,[1,LLL]);
for i=1:LLL
    modulated_noisy_picture(i)=modulated_noisy_picture(i)+noise(i);
end

demodulated_noisy_picture=demodulation(modulated_noisy_picture,fc,fs,Ts,M);

rcvd_noisy_symbols=detector(demodulated_noisy_picture,Sm,fc,fs,Ts,M);

rcvd_noisy_symbols=join(join(rcvd_noisy_symbols)'); rcvd_noisy_symbols=erase(rcvd_noisy_symbols," ");
rcvd_noisy_symbols=rcvd_noisy_symbols{1,1}(1:len_pic_str);

 
noisy_pic_array = decoder(rcvd_noisy_symbols,code,grayscales,size_a);
error_num=0;
for i=1:length(rcvd_noisy_symbols)
    if rcvd_noisy_symbols(i)~=rcvd_finaly_string(i)
        error_num=error_num+1;;
    end
end

%error_probability(cntr)=error_num/last_len;
%end
error_prob=error_num/last_len;
%figure
%plot(varia,error_probability*100,'blue')
%title('probability error based on variance');
%xlabel('variance');
%ylabel('error_probabilty(%)')

noisy_pic_array=uint8(noisy_pic_array);     %felan niazi besh nadarim
  nnoisy_pic_array=imresize(noisy_pic_array,8);    %felan niazi besh 
  %%
  figure
 imshow(nnoisy_pic_array);
 title('noisy 29.gif with variance=120');
%% Part E
[dem_r,~]=size(demodulated_picture); 
 corr_taw0_with_S0=zeros( 1 , dem_r );
corr_taw0_with_S1=zeros( 1 , dem_r );

 corr_taw0_with_S0_noisy=zeros( 1 , dem_r );
corr_taw0_with_S1_noisy=zeros( 1 , dem_r );

  ttemp_corr=zeros(1,2*Ts*fs+1);
 for i=1:dem_r
    ttemp_corr=xcorr(Sm(1,:),demodulated_picture(i,:));
    corr_taw0_with_S0(i)=ttemp_corr(Ts*fs+1);
    
    ttemp_corr=xcorr(Sm(2,:),demodulated_picture(i,:));
    corr_taw0_with_S1(i)=ttemp_corr(Ts*fs+1);
    
        ttemp_corr=xcorr(Sm(1,:),demodulated_noisy_picture(i,:));
    corr_taw0_with_S0_noisy(i)=ttemp_corr(Ts*fs+1);
    
    ttemp_corr=xcorr(Sm(2,:),demodulated_noisy_picture(i,:));
    corr_taw0_with_S1_noisy(i)=ttemp_corr(Ts*fs+1);
    
 end
 
%  figure
%  s=scatter(corr_taw0_with_S1,corr_taw0_with_S0);
% s.LineWidth = 0.6;
% s.MarkerEdgeColor = 'b';
% s.MarkerFaceColor = [0 0.5 0.5];
% hold on
% p1=plot([-150000:10000:150000],[-150000:10000:150000]) ;
% p2=plot([150000:-10000:-150000],[-150000:10000:150000]);
% hold off
% legend([p1 p2], 'margin1','margin2')
% title('(correlation of all symbols with S1 ,correlation of  all symbols with S0) for M=4 without noise in 38.gif')
%  

figure
 ss=scatter(corr_taw0_with_S1_noisy, corr_taw0_with_S0_noisy);
ss.LineWidth = 0.6;
ss.MarkerEdgeColor = 'b';
ss.MarkerFaceColor = [0 0.5 0.5];
hold on
p1=plot([-150000:10000:150000],[-150000:10000:150000]) ;
p2=plot([150000:-10000:-150000],[-150000:10000:150000]);
hold off
legend([p1 p2], 'margin1','margin2')
 title('(correlation of all symbols with S1 ,correlation of all symbols with S0) for M=4 with noise and variance=40 in 38.gif')
 
 