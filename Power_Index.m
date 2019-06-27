function [p] = Power_Index(W,n,T_i)

    p = zeros(n,n);  % power index (coln i = player i ; row i :power index of p_i)
    
    w_sum = zeros(n,1);
                                % calculate power index 
    for i = 1 : n               % for each node i update
        [a,b,v] = find(W(i,:)); % find non-zero element in W in row i,  
                                % index: (a,b), value vector: v
        nz = nnz(W(i,:));       % nr. of non zero element in W in row i
        
        for i_p = 1 : nz        % for each node i_p 
            for i_t = 1 : nz     % loop over each position 
                for t = 1 : T_i  % repeat T_i times
                    
                    P_per = randperm(nz); 
                    pos = find(P_per == i_p);
                    P_per([i_t pos]) = P_per([pos i_t]); % swap node position
                    
                    for i_x = 1 : i_t
                        
                        if w_sum(i) < 1/2
                            w_sum(i) = w_sum(i) + v(P_per(i_x));
                        end
                        
                        if w_sum(i) < 1/2 && (i_x == i_t)
                            w_sum(i) = 0; % set 0
                            break;
                        end
                        
                        if (w_sum(i) >= 1/2) && (i_x < i_t) % when weighted median unique 
                            
                            w_sum(i) = 0; % set 0
                            break;
                        elseif (w_sum(i) >= 1/2) && (i_x == i_t)
                            p(i,b(P_per(i_t))) = p(i,b(P_per(i_t))) + 1;
                            w_sum(i) = 0; % set 0
                            break;
                        end
                                           
                    end
                    w_sum(i) = 0; % set 0
                end                    
            end    
        end
        
        p(i,:) = p(i,:) / (T_i * nz);  % each element divided by T
        
%         for t = 1 : T_i % Approximate with T sample permutation(ex..ABC,ACB,BAC,BCA,CAB,CBA)           
%             P_per = randperm(nz); % random permutation with nr: 1 ~ nz            
%             for x_i = 1:nz 
%                 
%                 if w_sum(i) < 1/2  
%                     w_sum(i) = w_sum(i) + v(P_per(x_i));
%                     
%                 % elseif w_sum(i) == 1/2 % when weighted median not unique
%                     % p(P_per(t,x_i),i) = p(P_per(t,x_i),i) + 1;
%                 end
%                 
%                 if w_sum(i) >= 1/2 % when weighted median unique 
%                     p(i,b(P_per(x_i))) = p(i,b(P_per(x_i))) + 1;
%                     w_sum(i) = 0;
%                     break;
%                 end
%             end            
%         end
%         
%         p(i,:) = p(i,:) / T;  % each element divided by T
    end
    
    %P = sum(p)'; % sum of each column in p matrix, 
             % define power as  P(i) = p(1,i) + p(2,i) +...+ p(n,i) 
    
end