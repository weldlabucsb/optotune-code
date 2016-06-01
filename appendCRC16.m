function myans=appendCRC16(byte_array)

[low,high]=CRC16(byte_array);
myans=[byte_array, low, high];

end

function [low_byte, hi_byte]=CRC16(decinput)
poly16='8005';
poly16=[1 convert2arr(dec2bin(hex2dec(poly16),16))];

poly=poly16;
register=zeros(1,length(poly)-1);

for kk=1:length(decinput)
    bin=convert2arr(dec2bin(decinput(kk),8));
    bin=flip(bin);
    message=bin;
    if kk==length(decinput)
       message=[message zeros(1,length(poly)-1)]; 
    end
    
    while ~isempty(message)
        pop=register(1);
        register=[register(2:end) message(1)];
        message=message(2:end);
        if pop==1
           register=register+poly(2:end);
           for jj=1:length(register)
              register(jj)= mod(register(jj),2);
           end
        end
    end
end

register=flip(register);
str=convert2str(register);
byte=bin2dec(str);
hi_byte=bin2dec(str(1:8));
low_byte=bin2dec(str(9:16));

function mystr=convert2str(arr)
mystr='';
for ii=1:length(arr)
   mystr=[mystr num2str(arr(ii))];
end
end

function arr=convert2arr(str)
arr=[];

for ii=1:length(str)
   arr=[arr str2num(str(ii))]; 
end
end

end