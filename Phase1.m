function stuff = Phase1(P0,T0,RH0,m_air,P_in,P_ex)

% note on inputs
% p0 in psia
% t0 in fahrenheit
% RH0 in % (<1)
% m_air in lb/s
% P_in in inches of H20
% P_ex in inches of H20

T0S = ((T0-32)*5/9); % degrees C
T0K = T0S + 273.15;
P0S = P0*6.89475; % kPa
P_inS = P_in*0.249174;
P_exS = P_ex*0.249174;
m_airS = m_air*0.453592;

Psat = XSteam('psat_T',T0S);

% calculate x_wv
P_wv = RH0*Psat;
x_wv = (.622*P_wv/(P0-P_wv))/(1-(.622*P_wv/(P0-P_wv)));
y_wv = P_wv/P0;

% generate mole fractions
y_o2 = .150723*(1-y_wv);
y_n2 = .735258*(1-y_wv);
y_co2 = .053356*(1-y_wv);
y_ar = .012508*(1-y_wv);

whys = [y_o2,y_n2,y_co2,y_wv,y_ar]

% find intital properties
props1 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P0S,T0K);
props2 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, (P0S-P_inS),T0K);

T25S = binarysearch(whys,T0K,P0S,6)
