function mTemp = enthalpy_search(y,T,P,h_targ)

basis = prop_calc(y(1),y(2),y(3),y(4),y(5),P,T);
maxTemp = 4.*T;
minTemp = 0;    
% h = basis.h;
% P25 = P2*r_p;
%     h25s = basis.h;
%     %This assumes we know h1 already
%     h25a = (h25s-h1)/0.9 + h1;

while 1==1
    mTemp = (maxTemp+minTemp)/2;
    test = prop_calc(y(1),y(2),y(3),y(4),y(5),P,mTemp);
    if abs(h_targ-test.h) < .00000001
        break;
    end
    if h_targ > test.h
        minTemp = mTemp;
    end
    if h_targ < test.h
        maxTemp = mTemp;
    end
end
disp(mTemp);
result = mTemp;
