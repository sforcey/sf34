function [A,b,Aeq,beq,lb,ub] = BMEineq(n)
%This function determines the Matrix of vector inequalities that can be 
%used for LP and Genetic Programming Problem. LHS=A and RHS=b.
%lb is the lower bound on the x values for the problem while
%ub is the upper bound on the x values.

%Initialize Matrix to store constraints
k= nchoosek(n,2);
r = k+k*(n-2)+2*n; %+number of splits
A = zeros(r,k);
b = zeros(r,1);
Aeq = zeros(n,k); %What's the size of Aeq?
beq = zeros(n,1); %Size of beq?

%Create a vector for the lb on index:
lb = ones(k,1); 

%Create a vector for the ub on index:
ub = zeros(k,1);
for jj=1:k
    ub(jj) = 2^(n-3);
end

%Intersecting Cherry Facets: X_ik + X_jk - X_ij <= 2^(n-2)
%Start a counter to index the row
t = 1;
for i = 1:n-1
    for j = 1+i:n 
        for s = 1:n
            if i~=s && s~=j
                A(t,min(i,j)*(2*n-1-min(i,j))/2-n+max(i,j)) = -1;
                A(t,min(i,s)*(2*n-1-min(i,s))/2-n+max(i,s)) = 1;
                A(t,min(s,j)*(2*n-1-min(s,j))/2-n+max(s,j)) = 1;
                b(t,1) = 2^(n-3);
                t = t+1; %Increment t
            end  
        end
    end           
end

%Split Facets:sum(x_ij), where i,j are in S1 and i<j.
if n>5  %Condition for a split
    fl = floor(n/2);
    row = 0;
    for ii = 3:fl
        row = row+nchoosek(n,ii); 
    end
    col = fl;
    C = zeros(row,col); %Preallocates matrix for splits  
    x=1:n; %"Leaves" 
    for nn = 3:fl%Sizes of subsets
        %C = combnk(x,nn);%Generates a matrix of subsets for splits
        C1 = combnk(x,nn); %Generates a matrix of subsets for splits
        [rowC1,colC1] = size(C1);
        %Begin the transfer to C
        if nn == 3 
            C(1:rowC1,1:3) = C1(1:rowC1,1:3);
        else %nn>3 and we need to tack onto the end of the matrix
            C(all(~C,2),:) = [];
            [rC,cC] = size(C);  %Returns the current size of C
            C(rC+1:rC + rowC1,1:nn) = C1(1:rowC1,1:nn);
        end
    end
    
    %Generate the rows of A from the matrix of subsets
    d = 1; %Row counter 
    for rr = 1:row
        nonz = nonzeros(C(rr,:)); %Counts {|non-zero entries|}
        m = length(nonz);
        if mod(n,2)~=0 %n is odd, we want all size floor(n/2) subsets 
            %Enter as a split:
                for i = 1:n-1
                    for j=i+1:n
                        R1o = ismember(i,C(rr,:)); %Checks for i in C
                        R2o = ismember(j,C(rr,:)); %Checks for j in C
                        if R1o==1 && R2o==1 
                            A(k+k*(n-2)+2*n+d,i*(2*n-1-i)/2-n+j)=1;
                            b(k+k*(n-2)+2*n+d,1) = (m-1)*2^(n-3);
                        end
                    end
                end
                d =d+1; %Increment d for next row entry
         else %n is even so subsets will be slightly different     
               if m < n/2
                    for i = 1:n-1
                        for j=i+1:n
                            R1e = ismember(i,C(rr,:)); 
                            R2e = ismember(j,C(rr,:));
                            if R1e==1 && R2e==1 
                                A(k+k*(n-2)+2*n+d,i*(2*n-1-i)/2-n+j)=1;
                                b(k+k*(n-2)+2*n+d,1)=(m-1)*2^(n-3);
                            end
                        end
                    end
                    d =d+1; %Increment d for next entry
               else %m = n/2
                   for i = 1:n-1
                        for j=i+1:n
                            R1ef = ismember(1,C(rr,:));
                            R2ef = ismember(i,C(rr,:));
                            R3ef = ismember(j,C(rr,:));
                            if R1ef==1 && R2ef==1 && R3ef==1 
                                A(k+k*(n-2)+2*n+d,i*(2*n-1-i)/2-n+j)=1;
                                b(k+k*(n-2)+2*n+d,1) = (m-1)*2^(n-3);
                            end
                        end
                   end
                   d =d+1; %Increment d for next entry
                end    
        end
    end
end

%Equalities
%Dimension Restricting Equalities: sum(X_ij) = 2^(n-2)
%Since these are an equality,put them into Aeq and beq
for i = 1:n
    for j = 1:n 
        if i~=j 
        Aeq(i,min(i,j)*(2*n-1-min(i,j))/2-n+max(i,j)) = 1;
        end
    end
    beq(i,1) = 2^(n-2);
end
%Clean up and filter out zero rows in A and b
     A(all(~A,2),:) = [];
     b(all(~b,2),:) = [];
     Aeq(all(~Aeq,2),:) = [];
     beq(all(~beq,2),:) = [];
end