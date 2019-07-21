delta_t	= 0.25; % sampling time

A = [1     delta_t     0                   0;
      0     1           -delta_t*K_2        0;
      0     0           1                   delta_t;
      0     0           -delta_t*K_1*K_pp  1-delta_t*K_1*K_pd];
  
B = [0; 0; 0; delta_t*K_1*K_pp];

Q = diag([50 1 10 1]);

R = 1;


[K, S, e] = dlqr(A, B, Q, R);
