clear all
close all

n = 9;  %Number of taxa
[A,b,Aeq,beq,lb,ub] = BMEineq(n); %Generate BME(n)
k = nchoosek(n,2); %Number of decision variables
M=1:k; %Number of variables required to be discrete
e=1e-4; %Tolerance parameter
maxiteration = 15000; %Maximum number of iterations 
noise = 1; %Determines if we want to add noise

%% Cases
%Provide a true solution and distance data for experiments

%n = 6:
    %Test: N6T1
    
%     xtrue = [2 4 4 4 2 1 1 4 8 8 2 1 2 1 4];
%     
%     d = [4 3 3 3 4 5 5 3 2 2 4 5 4 5 3];
    
    %Test: N6T2
    
%     xtrue = [8 1 4 2 1 1 4 2 1 2 4 8 4 2 4];
%     
%     d = [2 14 3 13 14 14 3 13 14 13 3 2 12 13 3];

    
%n = 7:
    %Test: N7T1
    
%     xtrue = [1 1 2 16 4 8 16 8 1 4 2 8 1 4 2 2 8 4 4 8 8];
%     
%     d = [17 18 5 3 5 3 3 14 18 16 16 15 19 17 17 6 4 4 6 4 4];
    
    %Test: N7T2
    
%     xtrue = [4 8 8 4 4 4 2 2 4 4 16 16 2 2 2 2 2 2 16 4 4];
%     
%     d = [11 7 7 13 14 11 10 10 12 13 8 2 12 13 10 12 13 10 3 12 13];
   
    
%n = 8:
    %Test: N8T1
    
%     xtrue = [8 8 16 8 4 16 4 32 4 2 1 16 1 4 2 1 16 1 16 8 8 8 16 4 ...
%     16 2 32 2];
%     
%     d = [4 4 3 4 5 3 5 2 5 6 7 3 5 5 6 7 3 7 3 4 4 4 3 5 3 6 2 6];
    
    %Test: N8T2
    
%     xtrue = [32 16 4 4 4 4 4 16 4 4 4 2 2 8 8 8 4 4 32 8 4 4 8 4 4 ...
%     16 16 32];
%     
%     d = [5 9 15 16 15 18 19 8 14 15 14 17 18 10 11 10 13 14 9 12 15 16 ...
%     13 16 17 9 10 9];
    
    
%n = 9:
    %Test: N9T1
    
    xtrue = [16 8 32 64 4 2 1 1 32 32 16 16 8 4 4 16 8 32 16 8 8 32 8 4 2 2 ...
    4 2 1 1 32 16 16 32 32 64];

    d = [4 5 3 2 6 7 8 8 3 3 4 4 5 6 6 4 5 3 4 5 5 3 5 6 7 7 6 7 8 8 3 4 4 ... 
    3 3 2];

    %Test: N9T2
    
%     xtrue = [32 8 8 8 4 2 64 2 16 16 16 8 4 32 4 16 64 8 4 8 4 16 32 16 ...
%     8 16 8 4 8 4 32 4 32 2 64 2];
%     
%     d = [4 8 9 9 11 16 4 13 6 7 7 9 14 4 11 5 3 7 12 8 9 6 4 9 9 6 ...
%     8 13 9 10 9 11 6 16 5 13];


%n = 10:
    %Test: N10T1
    
%     xtrue = [128 64 32 16 8 4 2 1 1 64 32 16 8 4 2 1 1 64 32 16 8 4 2 2 ...
%     64 32 16 8 4 4 64 32 16 8 8 64 32 16 16 64 32 32 64 64 128];
% 
%     d = [2 3 4 5 6 7 8 9 9 3 4 5 6 7 8 9 9 3 4 5 6 7 8 8 3 4 5 6 7 7 ... 
%     3 4 5 6 6 3 4 5 5 3 4 4 3 3 2];

    %Test: N10T2
    
%     xtrue = [64 128 8 32 8 8 2 4 2 64 16 64 16 16 4 8 4 8 32 8 8 2 4 2 ...
%     32 128 32 8 16 8 32 32 8 16 8 32 8 16 8 32 64 32 64 128 64];
%     
%     d = [3 3 11 9 11 12 16 15 17 4 10 8 10 11 15 14 16 ...
%     12 10 12 13 17 16 18 6 2 5 9 8 10 6 7 11 10 12 5 9 8 10 8 7 9 3 3 4];
    

%% Solve the BME Problem and calculate topological distances
if noise == 0 %No noise added, so we have perfect data
    
    %Set aside a pool of 2 cores
    %parpool(2)

    spmd
    %Algorithm 2 (heuristic = 1)
    t2s = cputime;
    [x2,val2,status2]=DILP1(d,A,b,Aeq,beq,lb,ub,M,e,maxiteration,1);
    t2f = cputime - t2s;
    RF2 = RFmetric(x2,xtrue,n);
    
    %Algorithm 1 (heuristic = 0)
    t1s = cputime;
    [x1,val1,status1]=DILP1(d,A,b,Aeq,beq,lb,ub,M,e,maxiteration,0);
    t1f = cputime - t1s;
    RF1 = RFmetric(x1,xtrue,n);
    
    %l = distance(x,n); %if we want to draw the graph
    end
    
%% Solve the same problem, but add some noise. Do this if noise = 1
else %noise == 1

    %Create the distribution
    mu = 0; %Mean
    sigma = 1; %Standard Deviation
    h = 1; %Scales the perturbations
    
    delta = h*normrnd(mu,sigma,1,k); %Construct the perturbation vector
    
    L_inf = max(abs(delta)); %How much noise we add
    
    dpert = d + delta; %Perturb the objective function
    
    %Set aside a pool of 2 cores
    %parpool(2)
    
    spmd
    %Algorithm 2 (heuristic = 1)
    t2s = cputime;
    [x2pert,val2pert,pert2_status]=DILP1(dpert,A,b,Aeq,beq,lb,ub,M,e,maxiteration,1);
    t2f = cputime - t2s;
    RF2 = RFmetric(x2pert,xtrue,n);
    

    %Algorithm 1 (heuristic = 0)
    t1s = cputime;
    [x1pert,val1pert,pert1_status]=DILP1(dpert,A,b,Aeq,beq,lb,ub,M,e,maxiteration,0);
    t1f = cputime - t1s;
    RF1 = RFmetric(x1pert,xtrue,n);
    end
    %l_pert = distance(xpert,n); %if we want to draw the graph
end