
%% Initialization and model definition
init; % Change this to the init file corresponding to your helicopter

% sampling time
delta_t	= 0.25;

% Discrete time system model. x = [lambda r p p_dot e e_dot]'
A1 = [1    delta_t  0                   0                    0           0;
      0    1        -delta_t*K_2        0                    0           0;
      0    0        1                   delta_t              0           0;
      0    0        -delta_t*K_1*K_pp   1-delta_t*K_1*K_pd   0           0;
      0    0        0                   0                    1     delta_t;
      0    0        0           0   -delta_t*K_3*K_ep  1-delta_t*K_3*K_ed];
  
B1 = [0                 0; 
      0                 0;
      0                 0;
      delta_t*K_1*K_pp  0;
      0                 0;
      0                 delta_t*K_3*K_ep];
  
global N mx

  % Number of states and inputs
mx = size(A1,2); % Number of states (number of columns in A)
mu = size(B1,2); % Number of inputs(number of columns in B)
  



% Initial values
x1_0 = pi;                             % Lambda
x2_0 = 0;                              % r
x3_0 = 0;                              % p
x4_0 = 0;                              % p_dot
x5_0 = 0;                              % e
x6_0 = 0;                              % e_dot 
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]'; % Initial values

% Time horizon and initialization
N  = 40;                               % Time horizon for states
M  = N;                                % Time horizon for inputs
z  = zeros(N*mx+M*mu,1);               % Initialize z for the whole horizon
z0 = z;                                % Initial value for optimization

% Bounds
ul         = -Inf*ones(mu,1);          % Lower bound on states (no bound)
uu         = Inf*ones(mu,1);
ul(1) 	   = -30*pi/180;               % Lower bound on control
uu(1) 	   = 30*pi/180;                % Upper bound on control

xl      = -Inf*ones(mx,1);             % Lower bound on states (no bound)
xu      = Inf*ones(mx,1);              % Upper bound on states (no bound)
xl(3)   = ul(1);                       % Lower bound on state x3
xu(3)   = uu(1);                       % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = genBegr2(N,M,xl,xu,ul,uu);    % hint: gen_constraints
%vlb(N*mx+M*mu)  = 0;                  % We want the last input to be zero
%vub(N*mx+M*mu)  = 0;                  % We want the last input to be zero

Q1 = zeros(mx,mx);
Q1(1,1) = 1;                               % Weight on state x1

P1 = [1 0; 0 1];                           % Weight on input
Q = genQ2(Q1,P1,N,M);                      % Generate Q, hint: gen_q

% Generate c, this is the linear constant term in the QP
c = zeros((mx+mu)*N,1);            

%% Generate system matrixes for linear model
Aeq = genA2(A1,B1,N,mx,mu);               % Generate A, hint: gen_aeq
beq = zeros(mx*N,1);
beq(1:mx,1) = A1*x0;                       % Generate b



%% Solve QP problem with fmincon

z0(1) = x1_0;

%objective function

fun = @(Z) Z'*Q*Z;

% Non-linear contstraint
alpha = 0.2;
beta = 20;
lambda_t = 2*pi/3;

options = optimoptions('fmincon','Algorithm','sqp');
z0 = quadprog(Q,c,[],[],Aeq,beq, vlb, vub);
z = fmincon(fun,z0,[],[],Aeq,beq,vlb,vub,@nonlcon, options);

Q2 = diag([10 1 0 0 30 0]);
[K, ~, ~] = dlqr(A1, B1, Q2, P1);

%% Plotting
% Control input from solution
u1  = [z(N*mx+1:mu:N*mx+M*mu);z(N*mx+M*mu-1)]; 
u2  = [z(N*mx+2:mu:N*mx+M*mu);z(N*mx+M*mu)]; 

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

u1   = [zero_padding; u1; zero_padding];
u2   = [zero_padding; u2; zero_padding];
x1  = [pi*unit_padding; x1; zero_padding];
x2  = [zero_padding; x2; zero_padding];
x3  = [zero_padding; x3; zero_padding];
x4  = [zero_padding; x4; zero_padding];
x5  = [zero_padding; x5; unit_padding*z(N*mx-1)];
x6  = [zero_padding; x6; zero_padding];

t = (0:delta_t:delta_t*(length(x1)-1))';

optimal_u = [t, u1, u2];

optimal_x = [t, x1, x2, x3, x4, x5, x6];


figure(1)

subplot(211)
stairs(t,u1),grid
ylabel('u1')
title('Optimal input')
subplot(212)
stairs(t,u2),grid
ylabel('u2')
title('q = 1')


figure(2)

subplot(611)
plot(t,x1,'m',t,x1,'mo'),grid
ylabel('lambda')
title('Optimal trajectory')
subplot(612)
plot(t,x2,'m',t,x2','mo'),grid
ylabel('r')
subplot(613)
plot(t,x3,'m',t,x3,'mo'),grid
ylabel('p')
subplot(614)
plot(t,x4,'m',t,x4','mo'),grid
xlabel('tid (s)'),ylabel('pdot')
subplot(615)
plot(t,x5,'m',t,x5','mo'),grid
xlabel('tid (s)'),ylabel('e')
subplot(616)
plot(t,x6,'m',t,x6','mo'),grid
xlabel('tid (s)'),ylabel('edot')




