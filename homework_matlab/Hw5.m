%%  MAE 573 - Homework 5
%   Anthony Camarillo
%   Student ID: 013201459

clear
clc

%{
Question 1
Given [x1_dot x2_dot]' = [1 1; -2 3][x1 x2]' + [0 1]'*u
y = [2 2][x1 x2]'
change state variable to:
z1 = 3*x1 + 2*x2
z2 = 7*x1 + 5*x2
and write new equation for variables z1 and z2
%}
A1 = [1 1; -2 3];
B1 = [0 1]';
C1 = [2 2];
P1_inv = [3 2; 7 5];
P1 = inv(P1_inv);
A1_bar = P1_inv*A1*P1
B1_bar = P1_inv*B1
C1_bar = C1*P1
D1_bar = 0

%{
Question 2
Consider system x_dot = A*x + B*u, y = C*x
where
a) A = [0 2 0; 1 2 0; -1 0 1], B = [0 1 1]', C = [1 0 1]
b) A = [-6 1 0; -9 0 1; 0 0 0], B = [1 2 1]', C = [3 -2 1]
c) A = [-5 1/6 -1/6 0; -1/2 -16/3 1/3 1/2; -1/2 -1/3 -14/3 1/2; 0 -1/6 1/6 -5]
Find the transformation x = P*x_bar that transforms the state equations into DCF
or JCF
%}
A2_a = [0 2 0; 1 2 0; -1 0 1];
B2_a = [0 1 1]';
C2_a = [1 0 1];
[P2_a, A2_a_bar] = jordan(A2_a);
P2_a_inv = inv(P2_a);
P2_a
A2_a_bar
B2_a_bar = P2_a_inv*B2_a
C2_a_bar = C2_a*P2_a

A2_b = [-6 1 0; -9 0 1; 0 0 0];
B2_b = [1 2 1]';
C2_b = [3 -2 1];
[P2_b, A2_b_bar] = jordan(A2_b);
P2_b_inv = inv(P2_b);
P2_b
A2_b_bar
B2_b_bar = P2_b_inv*B2_b
C2_b_bar = C2_b*P2_b

A2_c = [-5 1/6 -1/6 0; -1/2 -16/3 1/3 1/2; -1/2 -1/3 -14/3 1/2; 0 -1/6 1/6 -5];
[P2_c, A2_c_bar] = jordan(A2_c);
P2_c
A2_c_bar

%{
Question 3
Transform below into CCF, show the transformation matrix P.
x_dot = [1 2 1; 0 1 3; 1 1 1]*x + [1 0 1]'*u
%}
syms s
A3 = [1 2 1; 0 1 3; 1 1 1];
B3 = [1 0 1]';
U3 = ctrb(A3,B3);
det(eye(3)*s-A3);
% From the above equation a2 = -3, a1 = -1, a0 = 3
M3 = [-1 -3 1; -3 1 0; 1 0 0];
P3 = U3*M3
P3_inv = inv(P3);
A3_bar = P3_inv*A3*P3
B3_bar = P3_inv*B3

%{
Question 4
x_dot = [2 0 0; 0 2 0; 0 3 1]*x + [0 1 0; 1 0 1]'*u
y = [1 0 0; 0 1 0]*x
Is the system completely controllable and completely observable?
%}
A4 = [2 0 0; 0 2 0; 0 3 1];
B4 = [0 1; 1 0; 0 1];
C4 = [1 0 0; 0 1 0];
U4 = ctrb(A4, B4);
V4 = obsv(A4, C4);
% rank of U4 is 3 and rank of V4 is 2, therefore the system is 
%completeley controllable but not completeley observable

%{
Question 5
Consider the system given by x_dot = [2 1 1; 5 3 6; -5 -1 -4]*x + [1 0 0]'*u
y = [1 1 2]*x
a) Determine which modes are controllable and/or observable using controllable
Kalman decomposition
b) Find the transfer function of the system and observe which modes appear as system
poles
%}
A5 = [2 1 1; 5 3 6; -5 -1 -4];
B5 = [1 0 0]';
C5 = [1 1 2];
U5 = ctrb(A5, B5);
V5 = obsv(A5, C5);
rank(U5);
rank(V5);
% From U's linearly independent columns we can find P_controllable as
P5_ctrb = [U5(:,1) U5(:,2) [0 0 1]'];
P5_ctrb_inv = inv(P5_ctrb);
% From V's linearly independent rows we can find P_observable_inverse as
P5_obsv_inv = [V5(1,:); V5(2,:); 0 0 1];
P5_obsv = inv(P5_obsv_inv);

A5_bar_ctrb = P5_ctrb_inv*A5*P5_ctrb;
B5_bar_ctrb = P5_ctrb_inv*B5;
C5_bar_ctrb = P5_ctrb;

A5_bar_obsv = P5_obsv_inv*A5*P5_obsv;
B5_bar_obsv = P5_obsv_inv*B5;
C5_bar_obsv = P5_obsv;
% Both x1_bar and x2_bar are observable and controllable and x3_bar is neither
% observable nor controllable

[num5, den5] = ss2tf(A5, B5, C5, 0);
Ys5 = zpk(num5, den5, 1)

%{
Question 6
Determine whether system is completely controllable and/or completely observable.
If not identify the uncontrollable or unobservable states
A = [2 1 0 0 0 0 0; 0 2 0 0 0 0 0; 0 0 -3 0 0 0 0; 0 0 0 -3 0 0 0;
0 0 0 0 4 1 0; 0 0 0 0 0 4 0; 0 0 0 0 0 0 2]
B = [1 0 1 1 0 1 -1; -1 1 0 1 0 0 1]'
C = [0 1 -1 -1 0 1 1; 0 0 -1 -1 1 0 -1]
%}

A6 = [2 1 0 0 0 0 0; 0 2 0 0 0 0 0; 0 0 -3 0 0 0 0; 0 0 0 -3 0 0 0;
0 0 0 0 4 1 0; 0 0 0 0 0 4 0; 0 0 0 0 0 0 2];
B6 = [1 0 1 1 0 1 -1; -1 1 0 1 0 0 1]';
C6 = [0 1 -1 -1 0 1 1; 0 0 -1 -1 1 0 -1];
U6 = ctrb(A6,B6);
V6 = obsv(A6, C6);
rank(U6);
rank(V6);
% From the rank of U and V, the system is fully controllable but not
% fully observable
% From observing the V matrix, the unobservable states are x9_dot, x11_dot
% and x13_dot
