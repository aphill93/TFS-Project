function stuff = Phase1(P0,T0,RH0,m_air,P_in,P_ex)

% note on inputs
% p0 in psia
% t0 in fahrenheit
% RH0 in % (<1)
% m_air in lb/s
% P_in in inches of H20
% P_ex in inches of H20

T0S = ((T0-32)*5/9); % degrees C
T0K = T0S + 273.15; % degrees K
P0S = P0*6.89475; % kPa
P_inS = P_in*0.249174; %kPa
P_exS = P_ex*0.249174; %kPa
m_airS = m_air*0.453592; % kg/s
eff_hpt = .9; % determined from base case
eff_lpt = .89; % determined from base case

Psat = XSteam('psat_T',T0S);

% calculate x_wv
P_wv = RH0*Psat;
x_wv = (.622*P_wv/(P0-P_wv))/(1-(.622*P_wv/(P0-P_wv)))
y_wv = P_wv/P0;

% generate mole fractions
y_o2 = .150723*(1-y_wv);
y_n2 = .735258*(1-y_wv);
y_co2 = .053356*(1-y_wv);
y_ar = .012508*(1-y_wv);

whys = [y_o2,y_n2,y_co2,y_wv,y_ar]'

% find intital properties
P0S
T0K
props1 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P0S,T0K)
P2S = (P0S-P_inS);
P2S
T0K
props2 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P2S,T0K)


P25S = P2S*6 % r_lpc = 6
T2K = T0K;
sprintf('HERE')
T25K_s = binarysearch(whys,T2K,P2S,P25S)
props25_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P25S, T25K_s);
h25S_a = (props25_s.h-props2.h)/0.82 + props2.h;
T25K_a = enthalpy_search(whys,T2K,P25S,h25S_a)
props25_a = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P25S, T25K_a)

P3S = P25S*4 % r_hpc = 4
T3K_s = binarysearch(whys,T25K_a,P25S,P3S);
props3_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P3S, T3K_s);
h3S_a = (props3_s.h-props25_a.h)/0.84 + props25_a.h;
T3K_a = enthalpy_search(whys,T25K_a,P3S,h3S_a)
props3_a = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P3S, T3K_a)

T4S = (2200-32)*(5/9); % degrees C
T4K = T4S + 273.15 % degrees K
P4S = P3S
props4 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P4S, T4K)

% base case- nominal 32MW, reverse engineer to find efficiences
% h48_a = props4.h - (props3_a.h - props2.h) % work_hpt = work_lpc + work_hpc
% props3_a.h - props2.h
% T48K = enthalpy_search(whys,T4K,P48S,h48_a)
% props48 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P48S, T48K)
% T48K_s = binarysearch(whys, T4K, P4S, P48S)
% props48_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P48S, T48K_s)
% eff_hpt = (props48_s.h - props4.h)/(props48.h - props4.h) % THIS IS BADDDD!!!!
% 
% h5S = (props48.h*1000 - (32000000/m_airS))/1000
% P5S = P0S + P_exS
% T5K = enthalpy_search(whys, T48K, P48S, h5S)
% props5 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P5S, T5K)
% T5K_s = binarysearch(whys, T48K, P48S, P5S)
% props5_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P5S, T5K_s)
% eff_lpt = (props48.h - props5_s.h)/(props48.h - props5.h) % 300% efficiency!

% actual case with turbine efficiencies from base case
P48S = 71*6.89475; % kPa
T48K_s = binarysearch(whys, T4K, P4S, P48S)
props48_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P48S, T48K_s)
h48S_a = props4.h + ((props48_s.h-props4.h)/eff_hpt);
T48K = enthalpy_search(whys, T4K, P4S, h48S_a)
props48_a = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P48S, T48K)

P5S = P0S + P_exS
T5K_s = binarysearch(whys, T48K, P48S, P5S)
props5_s = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P5S, T5K_s)
h5S_a = props48_a.h + ((props5_s.h-props48_a.h)/eff_hpt);
T5K = enthalpy_search(whys, T4K, P4S, h5S_a)
props5_a = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P5S, T5K)

T6K = T5K
P6S = P0S
props6 = prop_calc(y_o2, y_n2, y_co2, y_wv, y_ar, P6S, T6K)
