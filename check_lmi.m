clc;clear all;
%% Parameters
M = 10;m=5;l=1;g=9.8;mu=0.2;
k1 = mu/l; k2 = m*g;k3=-mu;k4=-(M+m)*g/l;k5=m*l;k6=-m;k7=-1/l;
x1_min = 0; x1_max = 2;
x2_min = -5; x2_max = 5;
x3_min = -pi/12; x3_max = pi/12;
x4_min = -pi/4 ; x4_max = pi/4;
pass_x1 = (x1_max-x1_min)/1000;
pass_x2 = (x2_max-x2_min)/1000;
pass_x3 = (x3_max-x3_min)/1000;
pass_x4 = (x4_max-x4_min)/1000;

%% Find min_max of premise variables
i = 1;
for x3=x3_min:pass_x3:x3_max
    z1(i) = 1/(M+m*sin(x3)^2);
    i = i + 1;
end
z1_max = max(z1); z1_min = min(z1);
z1s = [z1_max, z1_min];

i = 1;
for x3=x3_min:pass_x3:x3_max
    z2(i) = cos(x3);
    i = i + 1;
end
z2_max = max(z2); z2_min = min(z2);
z2s = [z2_max, z2_min];

i = 1;
for x3=x3_min:pass_x3:x3_max
    z3(i) = sin(x3)/x3;
    i = i + 1;
end
z3_max = max(z3); z3_min = min(z3);
z3s = [z3_max, z3_min];
i = 1;
for x3=x3_min:pass_x3:x3_max
    for x4=x4_min:pass_x4:x4_max    
        z4(i) = sin(x3)*x4;
        i = i + 1;
    end
end
z4_max = max(z4); z4_min = min(z4);
a_m = (z4_max*z1_max+z4_min*z1_min)/2;
a_r = (z4_max*z1_max-z4_min*z1_min)/2;
b_m = (z4_max*z1_max*z2_max+z4_min*z1_min*z2_min)/2;
b_r = (z4_max*z1_max*z2_max-z4_min*z1_min*z2_min)/2;

z1=0;z2=0;z3=0;z4=0;

for i=1:1:2
    z1 = z1s(i);
    for j =1:1:2
        z2 = z2s(j);
        for k=1:1:2
            z3 = z3s(k);
            A{k+2*(j-1)+2*2*(i-1)} = [0 1 0 0;
                                      0 k3*z1 k2*z1*z2*z3 k5*a_m;
                                      0 0 0 1;
                                      0 k1*z1*z2 k4*z3*z1 k6*b_m];
            B{k+2*(j-1)+2*2*(i-1)} = [0;z1;0;k7*z1*z2];
        end
    end
end
Ha = [0 1 0 0;0 0 0 1]'; Wa = [0 0 0 k5*a_r;0 0 0 k6*b_r];

%% LMI solver

nx = size(A{1},1); 
nu = size(B{1},2);    
Inu = eye(2);

P = sdpvar(nx,nx);
taura = sdpvar(2,2);
for i = 1:1:8
    N{i} = sdpvar(nu,nx,'full');
end

Con = [P>=0];

for i=1:1:8
    Gamma = blkvar;
    Gamma(1,1) = A{i}*P - B{i}*N{i} +(A{i}*P - B{i}*N{i})';
    Gamma(2,1) = Wa*P;
    Gamma(2,2) = -taura*Inu;
    Gamma = sdpvar(Gamma);
    Con = [Con, Gamma <=0];
end

for i=1:1:7
    for j=i+1:1:8
        Gammaij = blkvar;
        Gammaij(1,1) = A{i}*P-B{j}*N{i} + A{j}*P - B{i}*N{j} + (A{i}*P-B{j}*N{i} + A{j}*P - B{i}*N{j})';
        Gammaij(2,1) = 2*Wa*P;
        Gammaij(2,2) = -2*taura*Inu;
        Gammaij = sdpvar(Gammaij);
        Con = [Con, Gammaij <=0];
    end
end
%Check observability rank(obsv(A,C)) - rank(A)
options_sdp=sdpsettings;
options_sdp.solver='sedumi'; % sedumi sdpt3
% options_sdp.sdpt3.maxit  =50;
options_sdp.sedumi.maxiter  =50;
options_sdp.shift=1e-5; % next two lines: numerical parameters
options_sdp.verbose=1;  % =0: suppress screenoutput, =1: allow screen output

crit = 0;
solpb = solvesdp(Con,[],options_sdp);
[primal,dres] = checkset(Con);


if isequal(sum(primal>0),length(primal)) % see its significance in yalmiperror; checkset(pblmi)
  disp('The problem has been found FEASIBLE');


else
  disp('The problem has been found IIIIIIIINNNFEASIBLE');  
end


