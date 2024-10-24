% Exam 1
% Anthony Camarillo
% 013201459

% Question 1
%
% a) [x1_dot, x2_dot, x3_dot, x4_dot]' = 
% [[ x3 ],
%  [ x4 ],
%  [-(1/m*l)*cos(x1)*u-(g/l)*sin(x1)-((k*a^2)/m*l^2)*(sin(x1)-sin(x2))*cos(x1)],
%  [-(1/m*l)*cos(x2)*u-(g/l)*sin(x2)-((k*a^2)/m*l^2)*(sin(x2)-sin(x1))*cos(x1)]
% ]
syms x1 x2 x3 x4 u m g l k a
A1 = [ x3 ; x4 ;
  -(1/m*l)*cos(x1)*u-(g/l)*sin(x1)-((k*a^2)/m*l^2)*(sin(x1)-sin(x2))*cos(x1);
  -(1/m*l)*cos(x2)*u-(g/l)*sin(x2)-((k*a^2)/m*l^2)*(sin(x2)-sin(x1))*cos(x1)
    ]
%b
J1 = jacobian(A1, [x1,x2,x3,x4])
B1 = jacobian(A1, [u])
A1_linearize = subs(J1, [x1,x2,x3,x4], [0,0,0,0])
B1_linearize = subs(B1, [x1,x2,x3,x4], [0,0,0,0])

%c
% To determine stability, check eigenvalues of linearized matrix
eig_1 = eig(A1_linearize)
% From the above, two of the eigenvalues are positive and two of the
% eigenvalues are negative, therefore the system is not stable

% Question 2
% response y(t) is given by C*expm(A*t)*x(0)+C*integral(expm(A*t)*B*u(tau)dTau+D*u
syms t tau
A2 = [-1 1 0; 0 -1 0; 0 0 -2]
B2 = [0 1 1]'
C2 = [1 0 1]
x0 = [1 0 0]'
u2 = sin(t)
y2 = C2*expm(A2*t)*x0 + C2*int(expm(A2*(t-tau))*B2*u,tau,[0,t])

% Question 3
syms s
A3 = [-2 1 1; 0 -2 2; 0 0 -3]
B3 = [-1 -1 1]'
C3 = [1 2 4]
U3 = ctrb(A3, B3)
expand(det(s*eye(3)-A3))
% the above gives the equation s^3 + 7*s^2 + 16*s + 12
% so a0 = 12, a1 = 16, a2 = 7, a3 = 1
M3 = [16 7 1; 7 1 0; 1 0 0]
P3 = U3*M3
P3_inv = inv(P3)
A3_bar = P3_inv*A3*P3
B3_bar = P3_inv*B3
C3_bar = C3

% Question 4
A4 = [-1 1 0; 0 -1 0; 0 0 -2]
B4 = [0 1 1]'
C4 = [1 1 0]
U4 = ctrb(A4, B4)
V4 = obsv(A4, C4)
ctrb_rank_4 = rank(U4)
obsv_rank_4 = rank(V4)
% The controllability matrix for this system is full rank while the
% observability matrix for this system is only rank = 2. Therefore the
% system is completely controllable but not completely observable.

% Question 5
A5 = [2 2 2; 0 2 0; 0 0 -2]
B5 = [1 0; 0 0; 0 1]
C5 = [0 1 0; 0 0 1]

% a)
U5 = ctrb(A5, B5)
ctrb_rank_5 = rank(U5) % rank is 2 therefore 1 state is uncontrollable
% From the controllability matrix it can be seen that the first two linearly
% independent columns are [1 0 0]' and [0 0 1]', therefore to form a full rank
% matrix P we can use [0 1 0]'
P5_ctrb = [1 0 0; 0 0 1; 0 1 0]
P5_ctrb_inv = inv(P5_ctrb)
A5_ctrb_bar = P5_ctrb_inv*A5*P5_ctrb
B5_ctrb_bar = P5_ctrb_inv*B5
C5_ctrb_bar = C5*P5_ctrb
% it can be seen that x3_bar is the uncontrollable state corresponding
% to the eigenvalue 2
% b)
% Since the uncontrollable state is corresponding to the eigenvalue of 2,
% the system is not stabilizable; the uncontrollable eigenvalue must be negative
% c)
V5 = obsv(A5, C5)
obsv_rank_5 = rank(V5) % rank is 2 therefore 1 state is unobservable
% From the observability matrix it can be seen that the first two linearly
% independent rows are [0 1 0] and [0 0 1], therefore to form matrix P5_inv
% we can use [1 0 0]
P5_obsv_inv = [0 1 0; 0 0 1; 1 0 0]
P5_obsv = inv(P5_obsv_inv)
A5_obsv_bar = P5_obsv_inv*A5*P5_obsv
B5_obsv_bar = P5_obsv_inv*B5
C5_obsv_bar = C5*P5_obsv
% it can be seen that x3_bar is the unobservable state corresponding
% to the eigenvalue 2
% d)
% Since the unobservable state is corresponding to the eigenvalue of 2,
% the system is not detectable; the unobservable eigenvalue must be negative


