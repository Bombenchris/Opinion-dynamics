function [x_t,X,t_conv,P,P_final] = WM_Model(T,x_0,W,n)

w_sum = zeros(n,1);
X = zeros(n,1);
P = zeros(n,1);
P_final = zeros(n,1);
t_conv = 0;
x_st = x_0;

for t = 1:T   % loop over time step 
    t_conv = t;
    X = [X , x_0]; % keep track of change of x_t

    [x_s,I] = sort(x_0);     % I describes the arrangement of the elements of A into B along the sorted dimension
    i = randsample(n,1);    % individual i is picked uniformly at random

    if t >= 3000            % check convergence...
        x_d1 = x_0 - X(:,t-100);
        x_d2 = x_0 - X(:,t-200);
        x_d3 = x_0 - X(:,t-300);
        p = 1;
        if (norm(x_d1,p) == 0)&&(norm(x_d2,p) == 0)&&(norm(x_d3,p) == 0) 
            
            break;
        end         
    end

    for x_i = 1:n    % find the median weighted x(i)            
        if (w_sum(i) < 1/2) && (W(i,I(x_i)) ~= 0)    
            w_sum(i) = w_sum(i) + W(i,I(x_i));
        end
        
        if w_sum(i) == 1/2  % when weighted median not unique             
            
            while (w_sum(i) == 1/2)
                x_til = x_i - 1;  % check whether w = 0 inbetween. search x backward.
                w_sum(i) = w_sum(i) - W(i,I(x_til)); % until w_sum(i) < 1/2
            end
            
            if x_0(i) <= x_s(x_til) 
                x_0(i) = x_s(x_til); % update x_0(i) 
                P(I(x_i)) =  P(I(x_i)) + 1;
                 w_sum(i) = 0;
                break;
            elseif x_0(i) >= x_s(x_i)  
                x_0(i) = x_s(x_i); % update x_0(i)
                P(I(x_i)) =  P(I(x_i)) + 1;
                w_sum(i) = 0;
                break;
            else
                 % doesn't update x_0(i),  and w_ii = 0
                break;
            end
        end

        if w_sum(i) > 1/2
            if x_i > 1 
                x_0(i) = x_s(x_i); % update x_0(i)
                
                J = find(x_s == x_s(x_i));
                P(I(x_i)) =  P(I(x_i)) + 1;
                % P(I(J)) = P(I(J)) +1; % the same opinion all follows one times
                
                w_sum(i) = 0;
                break;
            elseif x_i == 1
                x_0(i) = x_s(x_i); % update x_0(i)  
                
                J = find(x_s == x_s(x_i));
                P(I(x_i)) =  P(I(x_i)) + 1;
                % P(I(J)) = P(I(J)) +1; % the same opinion all follows one times
                
                w_sum(i) = 0;
                break;
            end
        end
    end
end   

P = P / t_conv;
x_t = x_0;   

% second definition for social power. 
for i_p = 1 : n 
    count = sum(x_t == x_st(i_p));
    P_final(i_p) = count/n ; 
end

end
