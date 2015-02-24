%% MinSeg Project - EE 547 (PMP) - Winter 2015
% prepared by Paul Adams
%

%% Initialization
function minseg_project()
close all
digits(3);
format shortG
numerical_precision = 1e-9;
syms a x I_p m_p L r_w I_cmw m_w r_w T_m dx da dda ddx
syms k_t R V k_b 
g = 9.80665;    % [m/s^2]

%% Dynamical Model of the MinSeg Robot
%%
% <html> <h3> Equation of Motion </h3> </html>
E = [-(I_p + m_p*L^2), m_p*L*cos(a); ...
      m_p*L*r_w^2*cos(a), -(I_cmw + m_w*r_w^2 + m_p*r_w^2)];
x = [dda; ddx];
b = [T_m - m_p*L*g*sin(a); ...
     T_m*r_w + m_p*L*r_w^2*da^2*cos(a)];
render_latex([latex(E) latex(x) ' = ' latex(vpa(b, 5)) ], 12, 1.3)

%%
% <html> <h3> DC Motor Torque relation. </h3> </html>
T_m = k_t/R*V + (k_t*k_b)/(R*r_w)*dx + (k_t*k_b)/(R)*da;
render_latex(['T_m = ' latex(T_m)], 12, 1.3)

%%
% <html> <h3> State-space Matrices. </h3> </html>
A = [0, 1, 0, 0;
    (g*L*m_p*(I_cmw + (m_p + m_w)*r_w^2))/(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2), -(k_b*k_t*(I_cmw + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2)), 0, -(k_b*k_t*(I_cmw + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*r_w*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
    0, 0, 0, 1;
    (g*L^2*m_p^2*r_w^2)/(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2), -(k_b*k_t*r_w*(I_p + L*m_p*(L + r_w)))/(R*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2)), 0, -(k_b*k_t*(I_p + L*m_p*(L + r_w)))/(R*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2))];    
render_latex(['A = ' latex(vpa(A, 3))], 12, 1.3)

B = [0; 
    -(k_t*(I_cmw + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
    0;
    -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_cmw*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2))];
render_latex(['B = ' latex(vpa(B, 3))], 12, 1.3)
    
