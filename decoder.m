function final_pic_array = decoder(decoded_char,codes_string,grayscales,size_a)
R=length(decoded_char);
[R_code,~]=size(codes_string);
final_pic_array=zeros(size_a);
pic_row=1;pic_col=1;
number=1;
chemidoonam=1;
while (chemidoonam<64*64*R)

    if (length(decoded_char)~=0)&&(number<R) &&(strlength(codes_string(number))<length(decoded_char)) %&& (number<R_code)
       
         if decoded_char( 1:strlength(codes_string(number)) ) == convertStringsToChars(codes_string(number)) 
            final_pic_array(pic_row,pic_col)=grayscales(number);
            if length(decoded_char)>length(codes_string{number,1})
            decoded_char=decoded_char(length(codes_string{number,1})+1  :end);
            else
                decoded_char='';
            end
            pic_col=pic_col+1;
            if pic_col==65
                pic_row=pic_row+1;   pic_col=1;
            end
            number=1;
        else
           number=number+1;
        end
        
    else
        break;
    end
    chemidoonam=chemidoonam+1;
end

end