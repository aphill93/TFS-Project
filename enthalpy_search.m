ffunction result = binarysearch(y,T2,P2,pNext)

basis = prop_calc(y(1),y(2),y(3),y(4),y(5),P2,T2);
maxTemp = 4.*T2;
minTemp = 0;    
sNext = basis.Smix;
%P25 = P2*r_p;
%     h25s = basis.h;
%     %This assumes we know h1 already
%     h25a = (h25s-h1)/0.9 + h1;

while 1==1
    mTemp = (maxTemp+minTemp)/2;
 %  disp(mTemp)
    test = prop_calc(y(1),y(2),y(3),y(4),y(5),pNext,mTemp);
 %  disp(s25s)
 %  disp(test.Smix)
    if abs(sNext-test.Smix) < .00000001
        break;
    end
    if sNext > test.Smix
        minTemp = mTemp;
    end
    if sNext < test.Smix
        maxTemp = mTemp;
    end
end
%disp(mTemp);
result = mTemp;
