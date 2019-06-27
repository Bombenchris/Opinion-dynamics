function p = Power_generate(W,n)

    p = zeros(n,n);  
    syms x y
    
    t_ii = tic;
    for i = 1 : n
        
        [~,b,v] = find(W(i,:)); % for node i,, non zero weight
                                % index: (~,b), value vector: v
                                
        nz = nnz(W(i,:)); % nr. of non zero element in W in row i
        
        Z = [];
        Z = find(v > 0.5);
        if ~isempty(Z)   % weight lager than 0.5, 
            p(i,b(Z)) = 1; % set shapley shubik index = 1, break the loop
            continue;
        end
        
        f = 1;
        for p_i = 1:nz 
            f = f * (1 + x^(v(p_i)) * y) ;  % create generating function f 
        end
        
        t_i = tic;
        
        for x_i = 1:nz
            
%             f = 1; % initialize f(x,y)  
%             % power of x = weight , power of y = nr. of possible node
%             
%             for p_i = 1:nz 
%                 
%                 if p_i == x_i
%                     f = f*1;
%                     continue;
%                 end
%                 f = f * (1 + x^(v(p_i)) * y) ;  % create generating function f 
%             end
            
            f_p = f / (1 + x^(v(x_i)) * y);
                        % delect weight from node x_i
            f_p = sym2sympoly(f_p,{x,y});
            

            Pn = sympoly2polyn(f_p); % extract all information from the expanded f(x,y)
            
            [row,~] = find(Pn.ModelTerms(:,2) >= (0.5 - v(x_i)) & Pn.ModelTerms(:,2) < 0.5);
                        % find the sum of weight within range.
                        % some problem with sum of weight = 0.5 ,not unique
            sum = 0;
            
            
            for n_i = 1 : length(row)
                y_i = Pn.ModelTerms(row(n_i),3); % nr.of node
                c_i = Pn.Coefficients(row(n_i)); % coefficient
                sum = sum + c_i * factorial(y_i) * factorial(nz-y_i-1);
            end
                       
            p(i,b(x_i)) = sum/factorial(nz); 

        end
        toc(t_i) 
        disp(i)
        
    end
    toc(t_ii)
end