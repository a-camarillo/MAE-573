%% MAE 573 - Homework 4
% Anthony Camarillo
% Student ID: 013201459

clear
clc

%{
Question 1
Analyze the stability of the nonlinear system using Lyapunov's Linearization
Method
x_1dot = x_2
x_2dot = 2x_1 - 0.5(x_1)^3 - (1 + x_1)x_2
%}

% It can be seen that x_e = 0, for no input the linearized system
% can be given as delta_x_dot = A*delta_x
x_e = [0 0; -2 0; 2 0]';
eigen_values = [];
% iterate over all equilibrium points and find eigenvalues for each point
for i = 1:3
  v = x_e(:,i);
  A = [0 1; (2 - 1.5*(v(1))^2 - (v(2))), - 1 - v(2)];
  eigenvalues1 = [eigen_values eig(A)];
end
eigenvalues1
%%
% For the equilibrium point x_e = [0 0]', the eigenvalues are found to be
% 1 and -2, indicating this system is unstable
% 
% For the equilibrium points x_e = [2 0]' and x_e = [-2 0]', the eigenvalues
% are found to be -0.5 +/- 1.94i, indicating the system is asymptotically stable


%{
Question 2
Determine the stability of the following nonlinear systems about x_e = 0
using Lyapunov's Linearization Method
1) x_dot = -x^3
2) x_1dot = (x_1 - x_2)(x_1^2 + x_2^2 - 1), x_2dot = (x_1 + x_2)(x_1^2 + x_2^2 - 1)
%}

% 1) x_dot = -x^3
x_e = 0;
A2_1 = -3*(x_e)^2
% Since the above results in a zero value, the stability of the system
% cannot be determined

% 2) x_1dot = (x_1 - x_2)(x_1^2 + x_2^2 - 1), x_2dot = (x_1 + x_2)(x_1^2 + x_2^2 - 1)
syms x1 x2
A2_original = [(x1 - x2)*(x1^2 + x2^2 - 1); (x1 + x2)*(x1^2 + x2^2 - 1)];
J2 = jacobian(A2_original,[x1,x2]);
A2 = subs(J2, [x1, x2], [0, 0]);
eigenvalues2 = eig(A2)
% Since the eigenvalues of x_e = 0 result in lambda = -1 +/- i, the system is
% asymptotically stable

%{
Question 3
Given the pendulum equation without friction
theta_double_dot + sin*theta = u, where u(t) is the control torque input
a) Determine the state-space equations of the system
b) Linearize the system about [theta theta_dot]' = [0 0]' and express
it in the form x_dot = Ax + Bu, y = Cx + Du
c) Simulate the response with u = -[3 4]x, set initial conditions to
theta(0) = 1 and theta_dot(0) = 0
%}
% a) 
% State variables are given by x1 = theta, x2 = theta_double_dot
% Then x1_dot = x2 = theta_double_dot and x2_dot = u - sin(x1)
% x_dot = Ax + Bu where A = [x2;-sin(x1)], x = [theta theta_dot]'
% B = [0 1]', C = [1 0]  and D = [0]
% b)
A3_original = [x2; -sin(x1)];
J3 = jacobian([A3_original],[x1 x2]);
A3 = subs(J3, [x1, x2], [0, 0]);
eigenvalues3 = eig(A3)
% since the eigenvalues are found to be +/-i, the system is marginally stable
% at x_e = [0 0]'
% c)
% For u = -[3 4]x, the state-space representation is
% [x1_dot x2_dot]' = [x2 -sin(x1)-3(x1)-4(x1)]'
% subject to initial conditions theta(0)=x1(0)=1 and theta_dot(0)=x2(0)=0
% the system becomes [x1_dot x2_dot]' = [0 -3.8415]'
% For u = 0, the state-space representation is
% [x1_dot x2_dot]' = [x2 -sin(x1)]'
% subject to initial conditions theta(0)=x1(0)=1 and theta_dot(0)=x2(0)=0
% the system becomes [x1_dot x2_dot]' = [0 -0.8415]'

%{
Question 4
Consider the dynamic equation for an inverted pendulum
m*L^2*theta_double_dot+B*theta_dot-m*g*L*sin(theta) = tau(t)
where tau is the torque input and theta is the angular position.
Assume m,L,B=1 and g=9.81 m/s^2.
a) Determine the state equations
b) Linearize the system about [theta_e theta_e_dot]'=[0 0]'
and express it in the form x_dot=Ax+Bu, y=Cx+Du
%}
% a) Determine the state equations
% x1 = theta, x1_dot = x2 = theta_dot
% x2 = theta_dot, x2_dot = tau(t)-x2 - 9.81sin(x1)
% b) Linearize the system about x_e = [0 0]'
syms u
A4_original = [ x2; -sin(x1)-3*x1-4*x2+u ];
J4_A = jacobian([A4_original], [x1, x2])
J4_B = jacobian([A4_original], [u])
% delta_A is found as [0 1; -cos(x1)-3, -4]
% delta_B is found as [0 1]'
% about x_e = [0 0]'
% x_dot = [0 1; -4, -4]*delta_x + [0 1]'*delta_u

%{
Question 5
Euler equations of a rotating rigid spacecraft are given by:
J1*omega1_dot = (J2 - J3)*omega2*omega3 + u1
J2*omega2_dot = (J3 - J1)*omega3*omega1 + u2
J1*omega3_dot = (J1 - J2)*omega1*omega2 + u3
Torque inputs apply the feedback law u_i = -k_i*w_i where k_i are
positive constants. Verify local stability of the equilibrium point
(omega1,omega2,omega3) = (0,0,0)
%}
%% 
% Finding the state equations
% Dividing through by the respective Js the state equations can be given as
% w1_dot = ((J2-J3)/J1)*omega2*omega3 + u1/J1
% w2_dot = ((J3-J1)/J2)*omega3*omega1 + u2/J2
% w3_dot = ((J1-J2)/J3)*omega1*omega2 + u1/J3
% then finding the respective jacobians for Lyapunov's Linearization
% about (omega1, omega2, omega3) = (0,0,0)
% [omega1_dot, omega2_dot, omega3_dot]' = [-k1/J1 0 0; 0 -k2/J2 0; 0 0 -k3/J3]
% since the eigenvalues of the resulting matrix are the diagonal entries
% and both k and J are positive constants, all eigenvalues of the system
% about (0,0,0) will be negative meaning the system is locally stable

%{
Question 6
a) Obtain non-linear state space equations
b) Linearize system about (theta_e, theta_e_dot) = (0,0)
%}
% a)
% [x1_dot, x2_dot] = 
% [ x2; 
% ((-m*L*cos(x1)*u)/2(m+M) + (m*g*L/2)*sin(x1) - (m^2*L^2/(2*(m+M)))*(x2^2)*(sin(2*x1)))/((m*L^2/12) + ((m*L^2)/4)*(sin(x1)^2)+(((m*M*L^2)/4*(m+M))*(cos(x1)^2)))]
% b) Linearize system about (0,0)
syms x1 x2 m M L g u
A6_original = [ x2; ((-m*L*cos(x1)*u)/(2*(m+M)) + (m*g*L/2)*sin(x1) - (m^2*L^2/(2*(m+M)))*(x2^2)*(sin(2*x1)))/((m*L^2/12) + ((m*L^2)/4)*(sin(x1)^2)+(((m*M*L^2)/4*(m+M))*(cos(x1)^2)))];
J6 = jacobian([A6_original],[x1, x2]); 
A6 = subs(J6, [x1, x2], [0, 0])
