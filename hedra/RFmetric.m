function [total] = RFmetric(x1,x2,n)
%This function computes the Robinson-Foulds distance using characteristic
%functions to determine differences in the splits. It takes vector valued
%input as well as the number of taxa in the trees. Here, x1 and x2 
%represent two trees. 

%% 

if n < 4 
   fprintf = ('Error. The number of taxa must be at least 4.');
   return
end

%Initialize the sum for the RF metric

total = 0;

%% Generate the set of partitions

fl = floor(n/2);
    row = 0;
    for ii = 2:fl %Want sets of size 2 to floor(n/2)
        row = row + nchoosek(n,ii); %Calculates how many rows are needed 
    end
    col = fl; %Columns needed
    C = zeros(row,col); %Preallocates matrix for partitions   
    x = 1:n; %"Leaves" 
    for nn = 2:fl %Sizes of subsets
        
        C1 = combnk(x,nn); %Generates a matrix of subsets for splits
        [rowC1,~] = size(C1); %Find dimensions of C1
        
        %Begin the transfer to C
        if nn == 2 
            C(1:rowC1,1:2) = C1(1:rowC1,1:2);
        else %nn > 3 and we need to tack onto the end of the matrix
            C(all(~C,2),:) = [];
            [rC,cC] = size(C);  %Returns the current size of C
            C(rC+1:rC + rowC1,1:nn) = C1(1:rowC1,1:nn);
        end
    end
    
%% Check to see if splits are the same

%Sum over the splits

for rr = 1:row
    nonz = nonzeros(C(rr,:)); %Counts {|non-zero entries|}
    m = length(nonz);
    if mod(n,2)~=0 %n is odd, we want all size floor(n/2) subsets 
        %Enter as a split:
        %Initialize a local sum to change for each row of C
        sum1 = 0;
        sum2 = 0;  

            for i = 1:n-1
                for j=i+1:n

                    R1o = ismember(i,C(rr,:)); %Checks for i in C
                    R2o = ismember(j,C(rr,:)); %Checks for j in C

                    if R1o==1 && R2o==1 %It is in the split

                        %Update the sums
                        sum1 = sum1 + x1(i*(2*n-1-i)/2-n+j);
                        sum2 = sum2 + x2(i*(2*n-1-i)/2-n+j);

                    end
                end
            end

        %Now compute the characteristic function for this split

        if (sum1 == (m-1)*2^(n-3) && sum2 ~= (m-1)*2^(n-3)) || ...
           (sum2 == (m-1)*2^(n-3) && sum1 ~= (m-1)*2^(n-3))
       
            total = total + 1; %Characteristic is 1, so update the sum
        
        else
            
            total = total;
            
        end
       
     else %n is even so subsets will be slightly different
     %Enter as a split:
     %Initialize a local sum to change for each row of C
     sum1 = 0;
     sum2 = 0;  
    
           if m < n/2
                for i = 1:n-1
                    for j=i+1:n
                        R1e = ismember(i,C(rr,:)); 
                        R2e = ismember(j,C(rr,:));
                        if R1e==1 && R2e==1 
                            
                            %Update the sums
                            sum1 = sum1 + x1(i*(2*n-1-i)/2-n+j);
                            sum2 = sum2 + x2(i*(2*n-1-i)/2-n+j);
                            
                        end
                    end
                end

           else %m = n/2
               for i = 1:n-1
                    for j=i+1:n
                        R1ef = ismember(1,C(rr,:));
                        R2ef = ismember(i,C(rr,:));
                        R3ef = ismember(j,C(rr,:));
                        if R1ef==1 && R2ef==1 && R3ef==1
                            
                            %Update the sums
                            sum1 = sum1 + x1(i*(2*n-1-i)/2-n+j);
                            sum2 = sum2 + x2(i*(2*n-1-i)/2-n+j);
                            
                        end
                    end
               end
               
                       %Now compute the characteristic function for this split

               if (sum1 == (m-1)*2^(n-3) && sum2 ~= (m-1)*2^(n-3)) || ...
                    (sum2 == (m-1)*2^(n-3) && sum1 ~= (m-1)*2^(n-3))

                       total = total + 1; %Characteristic is 1, so update the sum

               else
            
                       total = total;                       
                       
               end
               
            end    
    end
end

end