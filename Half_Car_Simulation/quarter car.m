m1 = 290;
m2 = 59;
b1 = 1000;
k1 = 16182;
k2 = 19000;
A = [0 0 1 0;
    0 0 0 1;
    -k1/m1 k1/m1 -b1/m1 b1/m1;
    k1/m2 - (k1+k2)/m2 b1/m2 -b1/m2];
B = [0 0 0 k2/m2];
C = [1 0 0 0; 0 1 0 0];
D = [0 0]';
sys = ss(A, B, C, D);
t = linspace(0,100);
u = cos(t);
lsim(sys, u, t)