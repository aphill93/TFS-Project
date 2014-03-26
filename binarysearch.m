function result = binarysearch(y,T2,P2,r_p)

basis = prop_calc(y(1),y(2),y(3),y(4),y(5),P2,T2);
maxTemp = 2.*T2;
minTemp = 0;    
s25s = basis.Smix;
P25 = P2*r_p;
%     h25s = basis.h;
%     %This assumes we know h1 already
%     h25a = (h25s-h1)/0.9 + h1;

while 1==1
    mTemp = (maxTemp+minTemp)/2;
    test = prop_calc(y(1),y(2),y(3),y(4),y(5),P25,mTemp);
    if abs(s25s-test.Smix) < .00000001
        break;
    end
    if s25s > test.Smix
        minTemp = mTemp;
    end
    if s25s < test.Smix
        maxTemp = mTemp;
    end
end
result = mTemp;
