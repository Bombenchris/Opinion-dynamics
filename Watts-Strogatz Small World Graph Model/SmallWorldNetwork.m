%%
clc
clear
close all

n=1000; K=25; beta=0.2;

A = small_world(n, K, beta) + eye(n); 

%A = WattsStrogatz(n,K,beta) + eye(n);

A = A.*rand(n); % n-by-n matrix of random numbers, in the interval (0,1).
S = sum(A,2); % sum of the elements in each row of A
W = A.*(S.^-1); % row-stochastic influence matrix

% x_0 = 2.*rand(n,1) - 1; % initial opinion from uniform distribution [-1,1] at t = 0 ;
% x_0 = rand(n,1);% initial opinion from uniform distribution [0,1] at t = 0 ;
% x_0 = betarnd(2,2,n,1); % initial opinion obey symmetric uni-modal distribution given by a Beta distribution 
%% check nonzero weight 
nz = 0; 
for i = 1:n
    z = find(W(i,:));
    if size(z,2) >= nz
        nz = size(z,2);
    end
end


%%
W_inv = (W).^-1;
G_inv = digraph(W_inv);

W_inv_bet = (W+0.000001).^-1;
G_inv_bet = digraph(W_inv_bet);

%%
T = 100000; % Set update time step for WM-Model

% [x_t,X,t_conv] = WM_Model(T,x_0,W,n);  % WM_model dynamics

T_w = 1000; % about 600s for n = 1000 m = 2 
            % about 300s for n = 1000 m = 1 
P_W = zeros(n,1);
P_Final = zeros(n,1);

t = tic;
for i = 1:T_w
    x_0 = rand(n,1);  % uniform initial opinion from uniform distribution [0,1] at t = 0 ;
                        % need more evenly distributed for each loop.
                        % x_0(1) largest, second largest,...to
                        % smallest,each repeat T_e times 
    [x_t,X,t_conv,P_w,P_final] = WM_Model(T,x_0,W,n);
    P_W = P_W + P_w; 
    P_Final = P_Final + P_final; 
end
toc(t);

P_W = P_W/T_w;
P_Final = P_Final/T_w;

% [x_t,X,t_conv,P_W,P_final] = WM_Model(T,x_0,W,n);
%%

figure(1);
subplot(2,1,1)
histfit(x_0,20,'normal')
subplot(2,1,2)
histfit(x_t,20,'normal')

%subplot(3,1,3)
% plot(0:t_conv,X(100:150,1:t_conv+1));

%%
G = digraph(W); % graph to create an undirected graph or digraph to create a directed graph.
%% Generating function method (exact)
% p = Power_generate(W,n);  % shapley shubik power index network
% P_net = digraph(p);         % n = 1000, m=1 about s
% 
% P_inv = (P_net).^-1;
% P_inv_net = digraph(P_inv);

%% Structure sampling method (approximate)

T_p = 1000; % set time step for PI-Index calculation

t_1 = tic;                  % n = 1000, m=2 T_p = 1000 about 163s
p_I = Power_Index(W,n,T_p); % calculate power index
toc(t_1)

P_net_I = digraph(p_I);
%%
P_inv_I = (p_I).^-1;
P_inv_I_bet = (p_I+0.000001).^-1;

P_inv_I_net = digraph(P_inv_I);
P_inv_I_bet_net = digraph(P_inv_I_bet);
%%
c_deg = centrality(G,'indegree','Importance',G.Edges.Weight);
% define as c_deg(i) = sum over W(:,i)

c_pgrank = centrality(G,'pagerank','Importance',G.Edges.Weight);
% 'Importance' edge weights affect how the algorithm chooses successors. 
% Nodes with higher importance are more likely to be chosen

c_deg_p = centrality(P_net_I,'indegree','Importance',P_net_I.Edges.Weight); 
% define as c_deg(i) = sum over W(:,i)
c_pgrank_p = centrality(P_net_I,'pagerank','Importance',P_net_I.Edges.Weight);

% c_inclose = centrality(G_inv, 'incloseness','Cost', G_inv.Edges.Weight);
 [~, c_inclose, ~ ] = closenessCentrality( W_inv );
c_bet = centrality(G_inv_bet, 'betweenness','cost',G_inv_bet.Edges.Weight);

% c_inclose_p = centrality(P_inv_I_net, 'incloseness','Cost', P_inv_I_net.Edges.Weight);
 [~, c_inclose_p, ~ ] = closenessCentrality( P_inv_I );
c_bet_p = centrality(P_inv_I_bet_net, 'betweenness','cost',P_inv_I_bet_net.Edges.Weight);

%% for P_W
figure(2);

% compare with original network
b1 = c_deg\P_W;
b2 = c_pgrank\P_W;
b3 = c_deg_p\P_W;
b4 = c_pgrank_p\P_W;

ax1 = subplot(4,1,1);
r_1 = corrcoef(c_deg, P_W);
scatter(ax1,c_deg,P_W)
hold on
plot(c_deg,b1*c_deg)

xlabel(ax1,'Centrality: indegree(Orginal network)')
ylabel(ax1,'Power index')

% add correlation coeff. r 
str_1 = ['r_1= ',num2str(r_1(1,2)),';']; 
str_b1 = ['k_1= ',num2str(b1)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_1,str_b1)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

title('Social Power, define for each update step')

ax2 = subplot(4,1,2);
r_2 = corrcoef(c_pgrank, P_W);
scatter(ax2,c_pgrank,P_W)
hold on
plot(c_pgrank,b2*c_pgrank)

xlabel(ax2,'Centrality: pagerank(Orginal network)')
ylabel(ax2,'Power index')

% add correlation coeff. r 
str_2 = ['r_2= ',num2str(r_2(1,2)),';'];
str_b2 = ['k_2= ',num2str(b2)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_2,str_b2)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


% compare with shapley shubik network

ax3 = subplot(4,1,3);
r_3 = corrcoef(c_deg_p, P_W);
scatter(ax3,c_deg_p,P_W)
hold on
plot(c_deg_p,b3*c_deg_p)

xlabel(ax3,'Centrality: indegree(SSI Network)')
ylabel(ax3,'Power index')

% add correlation coeff. r 
str_3 = ['r_3= ',num2str(r_3(1,2)),';'];
str_b3 = ['k_3= ',num2str(b3)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_3,str_b3)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax4 = subplot(4,1,4);
r_4 = corrcoef(c_pgrank_p, P_W);
scatter(ax4,c_pgrank_p,P_W)
hold on
plot(c_pgrank_p,b4*c_pgrank_p)

xlabel(ax4,'Centrality: pagerank(SSI Network)')
ylabel(ax4,'Power index')

% add correlation coeff. r 
str_4 = ['r_4= ',num2str(r_4(1,2)),';'];
str_b4 = ['k_4= ',num2str(b4)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_4,str_b4)); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


%% for the second definition of social power P_Final

b1 = c_deg\P_Final;
b2 = c_pgrank\P_Final;
b3 = c_deg_p\P_Final;
b4 = c_pgrank_p\P_Final;

% original weighted network 
figure(3);

ax1 = subplot(4,1,1);
r_1 = corrcoef(c_deg, P_Final);
scatter(ax1,c_deg,P_Final)
hold on
plot(c_deg,b1*c_deg)

xlabel(ax1,'Centrality: indegree(original network)')
ylabel(ax1,'Power index')

% add correlation coeff. r 
str_1 = ['r_1= ',num2str(r_1(1,2)),';']; 
str_b1 = ['k_1= ',num2str(b1)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_1,str_b1)); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

title('Social Power, define for final step')

ax2 = subplot(4,1,2);
r_2 = corrcoef(c_pgrank, P_Final);
scatter(ax2,c_pgrank,P_Final)
hold on
plot(c_pgrank,b2*c_pgrank)

xlabel(ax2,'Centrality: pagerank(original network)')
ylabel(ax2,'Power index')

% add correlation coeff. r 
str_2 = ['r_2= ',num2str(r_2(1,2)),';'];
str_b2 = ['k_2= ',num2str(b2)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_2,str_b2)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

% compare in shapley shubik network 

ax3 = subplot(4,1,3);
r_3 = corrcoef(c_deg_p, P_Final);
scatter(ax3,c_deg_p,P_Final)
hold on
plot(c_deg_p,b3*c_deg_p)

xlabel(ax3,'Centrality: indegree(SSI network)')
ylabel(ax3,'Power index')

% add correlation coeff. r 
str_3 = ['r_3= ',num2str(r_3(1,2)),';'];
str_b3 = ['k_3= ',num2str(b3)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_3,str_b3)); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax4 = subplot(4,1,4);
r_4 = corrcoef(c_pgrank_p, P_Final);
scatter(ax4,c_pgrank_p,P_Final)
hold on
plot(c_pgrank_p,b4*c_pgrank_p)

xlabel(ax4,'Centrality: pagerank(SSI network)')
ylabel(ax4,'Power index')

% add correlation coeff. r 
str_4 = ['r_4= ',num2str(r_4(1,2)),';'];
str_b4 = ['k_4= ',num2str(b4)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_4,str_b4)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


%% for closeness and betweenness P_W

figure(4);
b1 = c_inclose'\P_W;
b2 = c_bet\P_W;
b3 = c_inclose_p'\P_W;
b4 = c_bet_p\P_W;

ax1 = subplot(4,1,1);
r_1 = corrcoef(c_inclose', P_W);
scatter(ax1,c_inclose',P_W)
hold on
plot(c_inclose',b1*c_inclose')

xlabel(ax1,'Centrality: inclose(Original Network)')
ylabel(ax1,'Power index')

% add correlation coeff. r 
str_1 = ['r_1= ',num2str(r_1(1,2)),';']; 
str_b1 = ['k_1= ',num2str(b1)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_1,str_b1)); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

title('Social Power, define for each update step')

ax2 = subplot(4,1,2);
r_2 = corrcoef(c_bet, P_W);
scatter(ax2,c_bet,P_W)
hold on
plot(c_bet,b2*c_bet)

xlabel(ax2,'Centrality: inbetw(Original Network)')
ylabel(ax2,'Power index')

% add correlation coeff. r 
str_2 = ['r_2= ',num2str(r_2(1,2)),';'];
str_b2 = ['k_2= ',num2str(b2)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_2,str_b2)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

ax3 = subplot(4,1,3);
r_3 = corrcoef(c_inclose_p', P_W);
scatter(ax3,c_inclose_p',P_W)
hold on
plot(c_inclose_p',b3*c_inclose_p')

xlabel(ax3,'Centrality: inclose_p(SSIN)')
ylabel(ax3,'Power index')

% add correlation coeff. r 
str_3 = ['r_3= ',num2str(r_3(1,2)),';'];
str_b3 = ['k_3= ',num2str(b3)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_3,str_b3)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax4 = subplot(4,1,4);
r_4 = corrcoef(c_bet_p, P_W);
scatter(ax4,c_bet_p,P_W)
hold on
plot(c_bet_p,b4*c_bet_p)

xlabel(ax4,'Centrality: inbetw_p(SSIN)')
ylabel(ax4,'Power index')

% add correlation coeff. r 
str_4 = ['r_4= ',num2str(r_4(1,2)),';'];
str_b4 = ['k_4= ',num2str(b4)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_4,str_b4)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

%% for closeness and betweenness P_Final

figure(5);

b1 = c_inclose'\P_Final;
b2 = c_bet\P_Final;
b3 = c_inclose_p'\P_Final;
b4 = c_bet_p\P_Final;

ax1 = subplot(4,1,1);
r_1 = corrcoef(c_inclose', P_Final);
scatter(ax1,c_inclose',P_Final)
hold on
plot(c_inclose',b1*c_inclose')

xlabel(ax1,'Centrality: inclose(original network)')
ylabel(ax1,'Power index')

% add correlation coeff. r 
str_1 = ['r_1= ',num2str(r_1(1,2)),';']; 
str_b1 = ['k_1= ',num2str(b1)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_1,str_b1)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

title('Social Power, define for final step')

ax2 = subplot(4,1,2);
r_2 = corrcoef(c_bet, P_Final);
scatter(ax2,c_bet,P_Final)
hold on
plot(c_bet,b2*c_bet)

xlabel(ax2,'Centrality: inbetw(original network)')
ylabel(ax2,'Power index')

% add correlation coeff. r 
str_2 = ['r_2= ',num2str(r_2(1,2)),';'];
str_b2 = ['k_2= ',num2str(b2)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_2,str_b2)); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

ax3 = subplot(4,1,3);
r_3 = corrcoef(c_inclose_p', P_Final);
scatter(ax3,c_inclose_p',P_Final)
hold on
plot(c_inclose_p',b3*c_inclose_p')

xlabel(ax3,'Centrality: inclose_p(SSIN)')
ylabel(ax3,'Power index')

% add correlation coeff. r 
str_3 = ['r_3= ',num2str(r_3(1,2)),';'];
str_b3 = ['k_3= ',num2str(b3)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_3,str_b3));  
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax4 = subplot(4,1,4);
r_4 = corrcoef(c_bet_p, P_Final);
scatter(ax4,c_bet_p,P_Final)
hold on
plot(c_bet_p,b4*c_bet_p)

xlabel(ax4,'Centrality: inbetw_p(SSIN)')
ylabel(ax4,'Power index')

% add correlation coeff. r 
str_4 = ['r_4= ',num2str(r_4(1,2)),';'];
str_b4 = ['k_4= ',num2str(b4)]; 
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),strcat(str_4,str_b4)); 

set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');
