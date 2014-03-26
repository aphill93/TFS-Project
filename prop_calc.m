function props = prop_calc(y_o2, y_n2, y_co2, y_h2o, y_ar, P, T)

% prop_calc computes thermodynamic properties  M, R, c_p, c_v, k, u,
% h, s, and p^o (all in mass based SI units) for a general mixture of gases
% containing O2, N2, CO2, H2O, and Argon.
%
% Inputs: mole fractions of each of the gases and the pressure (kPa) and
%         temperature (K) of the mixture
%
% Outputs: props- a struct containing M, R, c_p, c_v, k, u, h, s, and p^o

%y_o2, y_n2 = input('Please input each mole fraction followed by a comma:')
% sprintf('Oxygen: %d\nNitrogen: %d\nCarbon Dioxide: %d\nWater Vapor: %d\nArgon: %d'...
%         , y_o2, y_n2, y_co2, y_h2o, y_ar)

% inits and constants    
props = struct();
R_bar = 8.31441; % kpa*m^3/(kmol*K)
M_o2 = 31.999; R_o2 = R_bar/M_o2;
M_n2 = 28.013; R_n2 = R_bar/M_n2;
M_co2 = 44.01; R_co2 = R_bar/M_co2;
M_h2o = 18.015; R_h2o = R_bar/M_h2o;
M_ar = 39.948; R_ar = R_bar/M_ar;

% Reference values
T_ref = 298.15; P_ref = 101.325; % K and kPa
u_ref_o2 = 194.2; h_ref_o2 = 271.72; so_ref_o2 = 6.4107; % kJ/kg & kJ/kg*K
u_ref_n2 = 221.44; h_ref_n2 = 309.99; so_ref_n2 = 6.8405; % kJ/kg & kJ/kg*K
u_ref_co2 = 156.57; h_ref_co2 = 212.93; so_ref_co2 = 4.8585; % kJ/kg & kJ/kg*K
u_ref_h2o = 412.05; h_ref_h2o = 549.75; so_ref_h2o = 10.423; % kJ/kg & kJ/kg*K
u_ref_ar = 0; h_ref_ar = 0; T_ref_ar = 0; % no internal energy at abs zero 
so_ref_ar = 3.876; % kJ/kg*K


% Actual Calcs
M = (M_o2*y_o2 + M_n2*y_n2 + M_co2*y_co2 + M_h2o*y_h2o + M_ar*y_ar);

x_o2 = y_o2*M_o2/M;
x_n2 = y_n2*M_n2/M;
x_co2 = y_co2*M_co2/M;
x_h2o = y_h2o*M_h2o/M;
x_ar = y_ar*M_ar/M;

R = (R_o2*x_o2 + R_n2*x_n2 + R_co2*x_co2 + R_h2o*x_h2o + R_ar*x_ar);

a = 3; b= 4; c=7; d=8;
c_p_o2 = (25.48 + (1.52e-2)*T + (-.7155e-5)*T^2 + (1.312e-9)*T^3)/M_o2;
c_p_n2 = (28.9 + (-.1571e-2)*T + (0.8081e-5)*T^2 + (-2.873e-9)*T^3)/M_n2;
c_p_co2 = (22.26 + (5.981e-2)*T + (-3.501e-5)*T^2 + (7.469e-9)*T^3)/M_co2;
c_p_h2o = (32.24 + (0.1923e-2)*T + (1.055e-5)*T^2 + (-3.595e-9)*T^3)/M_h2o;
c_p_ar = 0.523; %const b/c it's a monotomic ideal gas

c_p = (c_p_o2*x_o2 + c_p_n2*x_n2 + c_p_co2*x_co2 + c_p_h2o*x_h2o + c_p_ar*x_ar);

c_v_o2 = c_p_o2 - R_o2;
c_v_n2 = c_p_n2 - R_n2;
c_v_co2 = c_p_co2 - R_co2;
c_v_h2o = c_p_h2o - R_h2o;
c_v_ar = c_p_ar - R_ar;


c_v = (c_v_o2*x_o2 + c_v_n2*x_n2 + c_v_co2*x_co2 + c_v_h2o*x_h2o + c_v_ar*x_ar);
c_v_bonus = c_p-R;

k = c_p/c_v;

fun=@(x) (25.48 - R_o2+ (1.52e-2)*x + (-.7155e-5)*x.^2 + (1.312e-9)*x.^3);
u_o2 = integral(fun,T_ref,T)/M_o2 + h_ref_o2;
fun2=@(x) (28.9 - R_n2 + (-.1571e-2)*x + (0.8081e-5)*x.^2 + (-2.873e-9)*x.^3);
u_n2 = integral(fun2,T_ref,T)/M_n2 + h_ref_n2;
fun3=@(x) (22.26 - R_co2 + (5.981e-2)*x + (-3.501e-5)*x.^2 + (7.469e-9)*x.^3);
u_co2 = integral(fun3,T_ref,T)/M_co2 + h_ref_co2;
fun4=@(x) (32.24 - R_h2o + (0.1923e-2)*x + (1.055e-5)*x.^2 + (-3.595e-9)*x.^3);
u_h2o = integral(fun4,T_ref,T)/M_h2o + h_ref_h2o;
u_ar = c_v_ar*(T - T_ref_ar) + u_ref_ar;

u = (u_o2*x_o2 + u_n2*x_n2 + u_co2*x_co2 + u_h2o*x_h2o + u_ar*x_ar);

fun=@(x) (25.48 + (1.52e-2)*x + (-.7155e-5)*x.^2 + (1.312e-9)*x.^3);
h_o2 = integral(fun,T_ref,T)/M_o2 + h_ref_o2;
fun2=@(x) (28.9 + (-.1571e-2)*x + (0.8081e-5)*x.^2 + (-2.873e-9)*x.^3);
h_n2 = integral(fun2,T_ref,T)/M_n2 + h_ref_n2;
fun3=@(x) (22.26 + (5.981e-2)*x + (-3.501e-5)*x.^2 + (7.469e-9)*x.^3);
h_co2 = integral(fun3,T_ref,T)/M_co2 + h_ref_co2;
fun4=@(x) (32.24 + (0.1923e-2)*x + (1.055e-5)*x.^2 + (-3.595e-9)*x.^3);
h_h2o = integral(fun4,T_ref,T)/M_h2o + h_ref_h2o;

h_ar = c_p_ar*(T - T_ref_ar) + h_ref_ar;

h = (h_o2*x_o2 + h_n2*x_n2 + h_co2*x_co2 + h_h2o*x_h2o + h_ar*x_ar);

so_o2 = c_p_o2*log(T/T_ref) + so_ref_o2;
so_n2 = c_p_n2*log(T/T_ref) + so_ref_n2;
so_co2 = c_p_co2*log(T/T_ref) + so_ref_co2;
so_h2o = c_p_h2o*log(T/T_ref) + so_ref_h2o;
so_ar = c_p_ar*log(T/T_ref) + so_ref_ar;

Smix1 = (so_o2*x_o2 + so_n2*x_n2 + so_co2*x_co2 + so_h2o*x_h2o + so_ar*x_ar);
Smix_ref = (so_ref_o2*x_o2 + so_ref_n2*x_n2 + so_ref_co2*x_co2 + so_ref_h2o*x_h2o + so_ref_ar*x_ar);

log_o2 = y_o2;
log_n2 = y_n2;
log_co2 = y_co2;
log_h2o = y_h2o;
log_ar = y_ar;

if(y_o2 == 0)
    log_o2 = 1;
end
if(y_n2 == 0)
    log_n2 = 1;
end
if(y_co2 == 0)
    log_co2 = 1;
end
if(y_h2o == 0)
    log_h2o = 1;
end
if(y_ar == 0)
    log_ar = 1;
end


Smix2 = (y_o2*R_o2*log(log_o2) + y_n2*R_n2*log(log_n2) + y_co2*R_co2*log(log_co2)...
    + y_h2o*R_h2o*log(log_h2o) + y_ar*R_ar*log(log_ar));

Smix3 = R*log(P/P_ref);

Smix = Smix1 - Smix2 - Smix3;


Po = (exp(Smix1/R)/exp(Smix_ref/R));
%Po_bonus = P/P_ref;

props = struct('M',M,'R',R,'c_p',c_p,'c_v',c_v,'k',k,'u',u,'h',h,'Smix',Smix,'Po',Po);


% fid = fopen('resultslast.xls','a')
% fprintf(fid, '%d\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n %f\n', T, M, R, c_p, c_v, k, u, h, Smix, Po );
% fclose(fid);
