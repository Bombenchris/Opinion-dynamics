%% Test centrality distance
n=3;
A = rand(n); % n-by-n matrix of random numbers, in the interval (0,1).
S = sum(A,2);
W = A.*(S.^-1);
G = digraph(W);

W_G = G.Edges.Weight;
% C_cl = centrality(G, 'incloseness','Cost',W_G);
C_cl = node_betweenness_slow(W);
h = plot(G,'MarkerSize',5);
h.NodeCData = C_cl;
colormap jet
colorbar


nrReachable = n - 1; % do not count starting node
distSum = sum(W(2:3,1));
c = 1 / distSum;
%% Test WM model
clc
clear
close all

W = [0.3, 0.3, 0.3,0.1; 
     0.3, 0.3, 0.3,0.1; 
     0.3, 0.3, 0.3,0.1;
     0.3, 0.3, 0.3,0.1]; % w_ii = 0.5, row stochastic
n= 4;

% W = [0.6, 0.2, 0.2; 
%      0.1, 0.7, 0.3; 
%      0.05, 0.15, 0.8]; % w_ii > 0.5, row stochastic
%  
% W = [0.3, 0, 0,0.55,0.15;
%     0.2, 0, 0.4,0.4,0;
%     0.3, 0.25, 0.25,0.2,0;
%     0.1, 0.15, 0.25,0.25,0.25;
%     0.6,0.1,0.1,0.1,0.1];
% n=10;
% A = rand(n); % n-by-n matrix of random numbers, in the interval (0,1).
% S = sum(A,2);
% W = A.*(S.^-1);

G = digraph(W);

p = Power_generate(W,n);  % shapley shubik power index network
P_net = digraph(p);

x_0 = rand(n,1);
T = 5000;

T_w = 100;
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

%[x_t,X,t_conv,P_W] = WM_Model(T,x_0,W,n);

%%
% test sampling model
% T_st = 5000;
% [P,p] = Power_Index(W,n,T_st);

c_deg = centrality(G,'indegree','Importance',G.Edges.Weight); 
% define as c_deg(i) = sum over W(:,i)
c_pgrank = centrality(G,'pagerank','Importance',G.Edges.Weight);

c_deg_p = centrality(P_net,'indegree','Importance',P_net.Edges.Weight); 
% define as c_deg(i) = sum over W(:,i)
c_pgrank_p = centrality(P_net,'pagerank','Importance',P_net.Edges.Weight);

%C_cl = centrality(G, 'incloseness','Cost', G.Edges.Weight);
%C_out = centrality(G, 'outdegree','Importance',G.Edges.Weight);
% outdegree wont make sense, cause sum(W(i,:)) = 1

ax1 = subplot(4,1,1);
r_1 = corrcoef(c_deg, P_W);
scatter(ax1,c_deg,P_W)
xlabel(ax1,'Centrality: indegree')
ylabel(ax1,'Power index')

% add correlation coeff. r 
str_1 = ['r_1= ',num2str(r_1(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),str_1); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax2 = subplot(4,1,2);
r_2 = corrcoef(c_pgrank, P_W);
scatter(ax2,c_pgrank,P_W)
xlabel(ax2,'Centrality: pagerank')
ylabel(ax2,'Power')

% add correlation coeff. r 
str_2 = ['r_2= ',num2str(r_2(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),str_2); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');

ax3 = subplot(4,1,3);
r_3 = corrcoef(c_deg_p, P_W);
scatter(ax3,c_deg_p,P_W)
xlabel(ax3,'Centrality: indegree')
ylabel(ax3,'Power index')

% add correlation coeff. r 
str_3 = ['r_3= ',num2str(r_3(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),str_3); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');


ax4 = subplot(4,1,4);
r_4 = corrcoef(c_pgrank_p, P_W);
scatter(ax4,c_pgrank_p,P_W)
xlabel(ax4,'Centrality: pagerank')
ylabel(ax4,'Power')

% add correlation coeff. r 
str_4 = ['r_4= ',num2str(r_4(1,2))];
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')),str_4); 
set(T, 'fontsize', 14, 'verticalalignment', 'top', 'horizontalalignment', 'left');



%% Test Power_generate function
clc
clear
close all

% W = [1, 0, 0; 
%      0.1, 0.9, 0; 
%      0.35, 0, 0.65]; % w_ii = 0.5, row stochastic
% n = 3 
% W = [0.6, 0.2, 0.2; 
%      0.1, 0.7, 0.3; 
%      0.05, 0.15, 0.8]; % w_ii > 0.5, row stochastic
%  
W = [1, 1, 1,0,0;
     1, 1, 0,0,0;
     1, 0, 1,1,1;
     1, 0, 1,0,0;
     4, 0, 1,0,0];

n = 5;

% W = [0.4, 0, 0,0;  % test incloseness
%      0.2, 0, 0,0;
%      0.2, 0, 0,0;
%      0,   0, 0.2,0];
% n = 4;

% W = [0.2, 0, 0;  % Test betweenness
%      0.2, 0, 0; 
%      0.2, 0, 0]; 
% 
% n = 3; 

% n=20;
% A = rand(n); % n-by-n matrix of random numbers, in the interval (0,1).
% S = sum(A,2);
% W = A.*(S.^-1);

% W = [0.25, 0.25, 0.25,0;  % test incloseness
%      0.25, 0.25, 0.25,0;
%      0.25, 0.25, 0.25,0;
%      0, 0.5, 0,0.5];
% 
% n = 4;

G = digraph(W);

% W = [0.3, 0.4, 0.3; 
%      0.6, 0.2, 0.2; 
%      0.2, 0.2, 0.6]; % w_ii = 0.5, row stochastic
% n = 3

% p = Power_generate(W,n);

%c_inclose = centrality(G_inv, 'incloseness','Cost', G_inv.Edges.Weight);
% C = closeness(W_inv);
 [~, inCloseness, ~ ] = closenessCentrality( W.^-1 );
%hc = harmoniccentrality(n,W_inv,'incloseness');

W_inv = (W).^-1;
G_inv = digraph(W_inv);

c_bet = centrality(G_inv,'betweenness','cost',G_inv.Edges.Weight);
% bc = betweenness_centrality(W_inv);
% betw = node_betweenness_faster(W_inv);
%%
T_p = 500; % set time step for PI-Index calculation
p_I = Power_Index(W,n,T_p); % calculate power index
