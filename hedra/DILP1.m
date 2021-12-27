%POLYSPLIT
%[x,val,status]=DILP1(f,A,b,Aeq,beq,lb,ub,M,e,maxiteration,heuristic)
%This function solves a discrete-integer linear programming problem
%using the branch and bound algorithm.
%The code uses MATLAB's linear programming solver "linprog"
%to solve the LP relaxations at each node of the branch and bound tree.
%   min f*x
%  subject to
%        A*x <=b
%        Aeq * x = beq   
%        lb <= x <= ub
%        M is a vector of indices for the discrete variables 
%        e is the tolerance
%The output variables are:
% x : the solution
% val: value of the objective function at the optimal solution
% status =1 if successful
%        =0 if maximum number of iterations reached in the linprog function
%        =-1 if there is no solution
%In order to run this code, you will need the contraints from BMEineq.m 
%and BNBtest.m.
%Author: William Sands
%This code was adapted from code provided by Kartik Sivaramakrishnan (see
%MILP.m) to solve a particular discrete-integer linear programming problem. 

function [x,val,status]=DILP1(f,A,b,Aeq,beq,lb,ub,M,e,maxiteration,heuristic)
%%
global count
global maxiter0
global maxiter1
h = heuristic; %This determines the rounding scheme for the branching
count=2;

if h == 0
     maxiter0 = maxiteration; %We don't round any of the variables
elseif h == 1
    maxiter0 = maxiteration; %We round only one variable at a time
    maxiter1= 2*maxiteration;
else
    fprintf('error: heuristic must equal 0 or 1.')
    return
end

options = optimset('display','off');
options.Algorithm = 'dual-simplex';
options.ConstraintTolerance = '1e-10';

bound=inf; %The initial bound is set to +ve infinity
%Solve the LP relaxation at the root node using MATLAB's linprog function
%Type "help linprog" for help with the linprog routine

%Solve the initial LP solution to determine feasiblility
[x0,val0,exitflag]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 

%If the LP is infeasible, then don't branch
if exitflag <= 0
    x = [];
    val = [];
    status = exitflag;
    return
end

if h == 1  %Determine number of cherries
    cherries = find(abs(x0(M)-ub(M)) <= e);
    newub = find(abs(x0(M)-ub(M)) > e);
    [row,~] = size(cherries);
    if row >=2 
        %Initialize Aeq and Beq.
            [rq,cq] = size(Aeq);
            Aeq = [Aeq;zeros(row,cq)]; %Same column size as A
            beq = [beq;zeros(row,1)];
            %Set the cherries as equalities in Aeq and beq
            for rr = 1:row
                Aeq(rq+rr,cherries(rr)) = 1;
                beq(rq+rr,1) = ub(cherries(rr));
            end
            %Reset the upper bounds since we found 'all' of the cherries
            ub(newub) = ub(cherries(1))-0.5*ub(cherries(1)); %Set new ub
            [x,val,status,b]=branch1(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound);
            
    else %Don't fix the cherries. Cherry cannot be forced.
        [x,val,status,b]=branch1(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound);
    end
    
else %Don't use the rounding heuristic
    [x,val,status,b]=branch0(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound);
end
end
%%
function [xx,val,status,bb]=branch0(f,A,b,Aeq,beq,lb,ub,x,v,M,e,bound)
global count
global maxiter0

options.Display = 'off';
options.Algorithm = 'dual-simplex';
options.ConstraintTolerance = '1e-8';

%Solve the LP relaxation at the current node
[x0,val0,status0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options);

%If the LP relaxation is infeasible,then PRUNE THE NODE BY INFEASIBILITY
%If the new objective value is worse, then prune by optimality
if status0<=0 || val0 > bound 
    %Return the input to linprog
    xx=x;
    val=v;
    status=status0;
    bb=bound;
    return;
end

%%
%If the solution to the LP relaxation is feasible in the DILP problem, then check the objective value of this 
%against the objective value of the best feasible/valid solution that has been obtained so far for the DILP problem.
%If the new feasible/valid solution has a lower objective value then update the bound
%Else PRUNE THE NODE BY OPTIMALITY

%Calculate tolerance and find branching variables.
[K,ind]= max(min(abs(x0(M)-2.^floor(log2(x0(M)))),abs(x0(M)-2.^ceil(log2(x0(M))))));

if K<e || count> maxiter0 %If we have a valid solution or exceed maxiter
    status=1;        
    if val0 < bound %The new feasible solution is an improvement 
        xx=x0; 
        val=val0;
        bb=val0;
    else
        xx=x;  %Return the input solution and most recent bounds
        val=v;
        bb=bound;
    end
    return
end

%%
%If we come here this means that the solution of the LP relaxation is not valid in the DILP problem.
%However, the objective value of the LP relaxation is lower than the current bound.
%So we branch on this node to create two subproblems.
%We will solve the two subproblems recursively by calling the same branching function.

%Select the branching variable. 
br_var=M(ind(1));  
br_value=x0(br_var); 
[~,c]=size(A);

%First LP problem with the added constraint that x_i <= floor(x_i),i=ind(1)
A1=[A ; zeros(1,c)];
A1(end,br_var)=1;
b1=[b;2^floor(log2(br_value))];

%Second LP problem with the added constraint that x_i >= ceil(x_i),i=ind(1)
A2=[A ;zeros(1,c)];
A2(end,br_var)=-1;
b2=[b; -2^ceil(log2(br_value))];

%%
%Solve the first LP problem
count = count+1; %+One for each subproblem being created
[x1,val1,status1,bound1]=branch0(f,A1,b1,Aeq,beq,lb,ub,x0,val0,M,e,bound);
status = status1;
if status1 >0 && bound1<bound %If the solution was successfull and gives a better bound
    xx=x1;
    val=val1;
    bound=bound1;
    bb=bound1;
else
    xx=x0;
    val=val0;
    bb=bound;
end

%Solve the second LP problem
count = count+1; %+One for each subproblem being created
[x2,val2,status2,bound2]=branch0(f,A2,b2,Aeq,beq,lb,ub,x0,val0,M,e,bound);
if status2 >0 && bound2<bound %If the solution was successful and gives a better bound
    status=status2;
    xx=x2;
    val=val2;
    bb=bound2;
end
end
%%
%
%Branching function using two routines. It firsts applies a typical branch
%and bound algorithm, setting equalities when a decision variable is within
%its tolerance of its nearest allowed discrete variable. Once the function
%reaches the maximum allowable iterations, it switches to searching for a 
%valid solution using the current solution and the available constraints.
%
%%
function [xx,val,status,bb]=branch1(f,A,b,Aeq,beq,lb,ub,x,v,M,e,bound)
global count
global maxiter0
global maxiter1

options.Display = 'off';
options.Algorithm = 'dual-simplex';
options.ConstraintTolerance = '1e-8';

%Solve the LP relaxation at the current node
[x0,val0,status0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options);

%If the LP relaxation is infeasible,then PRUNE THE NODE BY INFEASIBILITY
%If the new objective value is worse, then prune by optimality
if status0<=0 || val0 > bound 
    %Return the input to linprog
    xx=x;
    val=v;
    status=status0;
    bb=bound;
    return;
end

%%
%If the solution to the LP relaxation is feasible in the DILP problem, then check the objective value of this 
%against the objective value of the best feasible/valid solution that has been obtained so far for the DILP problem.
%If the new feasible/valid solution has a lower objective value then update the bound
%Else PRUNE THE NODE BY OPTIMALITY

%Calculate tolerance and find branching variables.
[K,ind]= max(min(abs(x0(M)-2.^floor(log2(x0(M)))),abs(x0(M)-2.^ceil(log2(x0(M))))));

if K<e || count> maxiter1 %If we have a valid solution or exceed maxiter
    status=1;        
    if val0 < bound %The new feasible solution is an improvement 
        xx=x0; 
        val=val0;
        bb=val0;
    else
        xx=x;  %Return the input solution and most recent bounds
        val=v;
        bb=bound;
    end
    return
end

%%
%If we come here this means that the solution of the LP relaxation is not valid in the DILP problem.
%However, the objective value of the LP relaxation is lower than the current bound.
%So we branch on this node to create two subproblems.
%We will solve the two subproblems recursively by calling the same branching function.

%Select the branching variable. 
br_var=M(ind(1));  
br_value=x0(br_var); 
[req,ceq]=size(Aeq);
[~,c]=size(A);

%Use the rounding heuristic (optional) to check for entries arbitrarily close to their
%discrete values. Don't set equalities for variables that are already powers of 2.

%We are using a heuristic
nonz = find(abs(x(M)-2.^round(log2(x0(M))))>0); %nonpower of 2 entries

if count > maxiter0  
   cRind = find(min(abs(x0(nonz)-2.^round(log2(x0(nonz)))))<e); %Smallest within tol
end

if count > maxiter0
    if ~isempty(cRind) %If a variable can be rounded
        Aeq = [Aeq; zeros(1,ceq)]; 
        beq = [beq;zeros(1)]; 
        Aeq(req+1,cRind(1)) = 1; 
        beq(req+1) = 2^round(log2(x0(cRind(1)))); %Round that variable
    end
end

%First LP problem with the added constraint that x_i <= floor(x_i),i=ind(1)
A1=[A ; zeros(1,c)];
A1(end,br_var)=1;
b1=[b;2^floor(log2(br_value))];

%Second LP problem with the added constraint that x_i >= ceil(x_i),i=ind(1)
A2=[A ;zeros(1,c)];
A2(end,br_var)=-1;
b2=[b; -2^ceil(log2(br_value))];

%%
%Solve the first LP problem
count = count+1; %+One for each subproblem being created
[x1,val1,status1,bound1]=branch1(f,A1,b1,Aeq,beq,lb,ub,x0,val0,M,e,bound);
status = status1;
if status1 >0 && bound1<bound %If the solution was successfull and gives a better bound
    xx=x1;
    val=val1;
    bound=bound1;
    bb=bound1;
else
    xx=x0;
    val=val0;
    bb=bound;
end

%Solve the second LP problem
count = count+1; %+One for each subproblem being created
[x2,val2,status2,bound2]=branch1(f,A2,b2,Aeq,beq,lb,ub,x0,val0,M,e,bound);
if status2 >0 && bound2<bound %If the solution was successful and gives a better bound
    status=status2;
    xx=x2;
    val=val2;
    bb=bound2;
end
end