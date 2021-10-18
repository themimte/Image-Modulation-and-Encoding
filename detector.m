function [rcvd_symbol] = detector(demodulated_symbols,Sm,fc,fs,Ts,M)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

[row,column]=size(demodulated_symbols);
[row_Sm,column_Sm]=size(Sm);

rcvd_symbol=strings([1 row]);

max_sm=zeros(row,1);
num_of_Sm=0;
dot_Sm_symbols=zeros(row_Sm,1);
num_bit=log2(M);

for m=1:row
    for k=1:row_Sm
        dot_Sm_symbols(k,1)=dot(demodulated_symbols(m,:),Sm(k,:));        
    end
    
     max_dot=find(dot_Sm_symbols==max(dot_Sm_symbols));
     num_of_Sm=max_dot-1;
     max_sm(m)=num_of_Sm;
     temp=de2bi(num_of_Sm); temp=flip(temp); temp=string(temp); temp=join(temp); temp=erase(temp," ");
     if strlength(temp) < num_bit
       for num = strlength(temp)+1:num_bit
           temp=strcat("0",temp);
       end
     end

  rcvd_symbol{1,m}=temp{1,1};  
end

end