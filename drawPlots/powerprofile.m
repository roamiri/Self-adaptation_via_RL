%%
clear;
clc;
%% rho = defined, IL and CL original
T1 = 'oct4/T1/pro_x1_IL_%d_%d.mat';
T2 = 'oct4/T1/pro_x2_IL_%d_%d.mat';
T3 = 'oct4/T1/pro_x1_CL_%d_%d.mat';
T4 = 'oct4/T1/pro_x2_CL_%d_%d.mat';
T5 = 'oct4/T1/pro_greedy_%d_%d.mat';
% [P_min_FUE_2, P_sum_FUE_2, P_mean_FUE_2, P_max_FUE_2 ] = powerProfile(T3);
%% rho = 1, IL and CL original
T1 = 'oct8/T2/pro_x1_IL_%d_%d.mat';
T2 = 'oct8/T2/pro_x2_IL_%d_%d.mat';
T3 = 'oct8/T2/pro_x1_CL_%d_%d.mat';
T4 = 'oct8/T2/pro_x2_CL_%d_%d.mat';
T5 = 'oct4/T1/pro_greedy_%d_%d.mat';
%% rho = 1, IL and CL reverse reward functions
T1 = 'oct8/T3/pro_x1_IL_%d_%d.mat';
T2 = 'oct8/T3/pro_x2_IL_%d_%d.mat';
T3 = 'oct8/T3/pro_x1_CL_%d_%d.mat';
T4 = 'oct8/T3/pro_x2_CL_%d_%d.mat';
T5 = 'oct4/T1/pro_greedy_%d_%d.mat';
%% rho=1, IL vs CL with ()^2 Reward function
T1 = 'oct8/T2/pro_x1_IL_%d_%d.mat';
T2 = 'oct8/T2/pro_x2_IL_%d_%d.mat';
T3 = 'oct8/T3/pro_x1_CL_%d_%d.mat';
T4 = 'oct8/T3/pro_x2_CL_%d_%d.mat';
T5 = 'oct4/T1/pro_greedy_%d_%d.mat';
%% rho=1, IL vs CL with ()^3 Reward function
T1 = 'oct8/T3/pro_x1_IL_%d_%d.mat';
T2 = 'oct8/T3/pro_x2_IL_%d_%d.mat';
T3 = 'oct8/T2/pro_x1_CL_%d_%d.mat';
T4 = 'oct8/T2/pro_x2_CL_%d_%d.mat';
T5 = 'oct4/T1/pro_greedy_%d_%d.mat';
%% rho=1, IL , ()^3 Reward function, different values for \gamma_k
T1 = 'oct9/T1/pro_x1_IL_%d_%d.mat'; %\gamma=0
T2 = 'oct8/T3/pro_x1_IL_%d_%d.mat'; %\gamma = 0.5
T3 = 'oct9/T2/pro_x1_IL_%d_%d.mat'; %\gamma= 1.0
%% IL RF=() vs RF=()^3 
T1 = 'oct8/T3/pro_x1_IL_%d_%d.mat'; %()^3
T2 = 'oct9/T3/pro_x1_IL_%d_%d.mat'; %()^1
T5 = 'oct4/T1/pro_greedy_%d_%d.mat'; %greedy
%% IL, ()^2 vs ()^3
T1 = 'oct8/T2/pro_x1_IL_%d_%d.mat'; %()^2
T2 = 'oct8/T2/pro_x2_IL_%d_%d.mat'; %()^2
T3 = 'oct8/T3/pro_x1_IL_%d_%d.mat'; %()^3
T4 = 'oct8/T3/pro_x2_IL_%d_%d.mat'; %()^3
T5 = 'oct4/T1/pro_greedy_%d_%d.mat'; %greedy
%% IL, exp()^2 vs ()^3
T1 = 'oct11/T1/pro_x1_IL_%d_%d.mat'; %()^exp
T2 = 'oct11/T1/pro_x2_IL_%d_%d.mat'; %()^exp
T3 = 'oct8/T3/pro_x1_IL_%d_%d.mat'; %()^3
T4 = 'oct8/T3/pro_x2_IL_%d_%d.mat'; %()^3
T5 = 'oct4/T1/pro_greedy_%d_%d.mat'; %greedy
%% IL, ICC vs ()^3
T1 = 'oct12/T1/pro_x1_IL_%d_%d.mat'; %ICC
T2 = 'oct12/T1/pro_x1_CL_%d_%d.mat'; %ICC
T3 = 'oct8/T3/pro_x1_IL_%d_%d.mat'; %()^3
T4 = 'oct8/T3/pro_x1_CL_%d_%d.mat'; %()^3
T5 = 'oct4/T1/pro_greedy_%d_%d.mat'; %greedy
%%
T0 = 'Sep2/T2/pro_x2_IL_%d_%d.mat';
T1 = 'Sep2/T2/pro_x3_IL_%d_%d.mat';
T2 = 'Aug26/T1/pro_x2_IL_%d_%d.mat';
T3 = 'Aug26/T1/pro_x3_IL_%d_%d.mat';
T4 = 'Aug23/T2/pro_greedy_%d_%d.mat';
%%
[P_min_FUE, P_sum_FUE, P_mean_FUE, P_max_FUE ] = powerProfile(T1);
[P_min_FUE_1, P_sum_FUE_1, P_mean_FUE_1, P_max_FUE_1 ] = powerProfile(T2);
[P_min_FUE_2, P_sum_FUE_2, P_mean_FUE_2, P_max_FUE_2 ] = powerProfile(T3);
[P_min_FUE_3, P_sum_FUE_3, P_mean_FUE_3, P_max_FUE_3 ] = powerProfile(T4);
[P_min_FUE_4, P_sum_FUE_4, P_mean_FUE_4, P_max_FUE_4 ] = powerProfile(T5);
%%
figure;
hold on;
grid on;
box on;
% plot(ones(1,40)*1.0, '--k', 'LineWidth',1);
plot(P_mean_FUE, '--ok', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
plot(P_min_FUE, '--ob', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(P_max_FUE, '--or', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 15]);
legend({'IL+X_1'},'Interpreter','latex','FontSize',12);
%%
figure;
hold on;
grid on;
box on;
% plot(ones(1,40)*1.0, '--k', 'LineWidth',1);
plot(P_mean_FUE_1, '--ok', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
plot(P_min_FUE_1, '--ob', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(P_max_FUE_1, '--or', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 15]);
legend({'IL+X_2'},'Interpreter','latex','FontSize',12);
%%
figure;
hold on;
grid on;
box on;
% plot(ones(1,40)*1.0, '--k', 'LineWidth',1);
plot(P_mean_FUE_2, '--dk', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
plot(P_min_FUE_2, '--db', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(P_max_FUE_2, '--dr', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 15]);
legend({'CL+X_1'},'Interpreter','latex','FontSize',12);
%%
figure;
hold on;
grid on;
box on;
% plot(ones(1,40)*1.0, '--k', 'LineWidth',1);
plot(P_mean_FUE_3, '--dk', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
plot(P_min_FUE_3, '--db', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
plot(P_max_FUE_3, '--dr', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 15]);
legend({'CL+X_2'},'Interpreter','latex','FontSize',12);
%%
figure;
hold on;
grid on;
box on;
errorbar(1:10, P_mean_FUE, P_max_FUE-P_min_FUE, '--or', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
errorbar(1:10, P_mean_FUE_1, P_max_FUE_1-P_min_FUE_1, '--ob', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
errorbar(1:10, P_mean_FUE_2, P_max_FUE_2-P_min_FUE_2, '--og', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
errorbar(1:10, P_mean_FUE_3, P_max_FUE_3-P_min_FUE_3, '--ok', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
ylim([0 25]);
%%
%%
figure;

subplot(4,1,1);
hold on;
grid on;
box on;
% plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, P_mean_FUE, P_max_FUE-P_min_FUE, '--or', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','r', 'MarkerEdgeColor','b');
 xlim([1 10]);
ylim([0 20.0]);
legend({'IL+$\mathcal{X}_1$'}, 'Interpreter','latex','FontSize',12);
% legend({'qos', 'old pathloss'},'Interpreter','latex','FontSize',12);
% legend({'qos','IL+$\mathcal{X}_2$'},'FontSize',10, 'Interpreter','latex');


subplot(4,1,2);
hold on;
grid on;
box on;
% plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, P_mean_FUE_1, P_max_FUE_1-P_min_FUE_1, '--ob', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','b', 'MarkerEdgeColor','b');
 xlim([1 10]);
ylim([0 20.0]);
% legend({'qos','proximity RF'},'FontSize',12);
legend({'$\mathcal{X}_2$'},'Interpreter','latex','FontSize',12);
% legend({'qos','CL+$\mathcal{X}_2$'},'FontSize',10, 'Interpreter','latex');

subplot(4,1,3);
hold on;
grid on;
box on;
% plot(1:10, ones(1,10)*.50, '--k', 'LineWidth',1);
errorbar(1:10, P_mean_FUE_2, P_max_FUE_2-P_min_FUE_2, '--og', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 20.0]);
% legend('X_3');
legend({ '$\mathcal{X}_3$'},'Interpreter','latex','FontSize',12);
% legend({'qos','CL+$\mathcal{X}_1$'},'FontSize',10, 'Interpreter','latex');

subplot(4,1,4);
hold on;
grid on;
box on;
% plot(1:10, ones(1,10)*0.50, '--k', 'LineWidth',1);
errorbar(1:10, P_mean_FUE_3, P_max_FUE_3-P_min_FUE_3, '--ok', 'LineWidth',1.3,'MarkerSize',2, 'MarkerFaceColor','k', 'MarkerEdgeColor','b');
xlim([1 10]);
ylim([0 20.0]);
legend({'$\mathcal{X}_4$'},'FontSize',10, 'Interpreter','latex');

supertitle('','FontSize',14, 'FontWeight','bold');
%%
% figure;
hold on;
grid on;
box on;
% plot( ones(1,16)*2.0, '--k', 'LineWidth',1 );
plot(P_sum_FUE, '--or', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','r');%, 'MarkerEdgeColor','b');
plot(P_sum_FUE_1, '--ob', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','b');%, 'MarkerEdgeColor','b');
plot(P_sum_FUE_2, '--dr', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r');%, 'MarkerEdgeColor','b');
plot(P_sum_FUE_3, '--db', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b');%, 'MarkerEdgeColor','b');
plot(P_sum_FUE_4, '--*k', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k');
% title('SUM capacity of FUEs','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Sum power (mWatt)','FontSize',12);%, 'FontWeight','bold');
xlim([1 10]);
ylim([4 350]);
% legend({'CL+$\mathcal{X}_1,\rho=1$','CL+$\mathcal{X}_2,\rho=1$', 'CL+$\mathcal{X}_1$', 'CL+$\mathcal{X}_2$', 'greedy'},'Interpreter','latex','FontSize',12);
% legend({'IL','CL', '$\rho$'},'Interpreter','latex','FontSize',12);
% legend({'\alpha_1','\alpha_2', '\alpha_3'},'FontSize',14, 'FontWeight','bold');
% legend({'X_1','X_2', 'X_3', 'greedy'},'FontSize',12);%, 'FontWeight','bold');
%%
%%
function [P_min_FUE, P_sum_FUE, P_mean_FUE, P_max_FUE ] = powerProfile(T)
    P_C = [];    
    P_min_FUE = [];
    P_sum_FUE = [];
    P_mean_FUE = [];
    P_max_FUE = [];
    P_FUE_Mat = cell(1,10);
    P_FUE_Mat_W = cell(1,10);
    p_fue = [];
    for i=1:10
        fprintf('FBS num = %d\t', i);
        maxmue = 0.;
        maxfue = 0.;
        mue_C = 0.;
        minfue = 0.;
        sumfue = 0.;
        p_fue_vec = zeros(1,i);
        Cnt = 0;
        lowCnt = 0;

        for j=1:500
    %         s = sprintf('Aug16/IL/pro_IL_77_%d_%d.mat',i,j);
            s = sprintf(T,i,j);
            filename = strcat(s);
            if exist(s)
                load(filename);
                FBS = QFinal.FBS;
                for kk=1:size(FBS,2)
                     p_fue(1,kk) = FBS{1,kk}.P;
                end
                p_fue_vec = p_fue_vec + p_fue;
                Cnt = Cnt+1;
            end
        end
        fprintf('Total Cnt = %d\n',Cnt);

        P_FUE_Mat{i} = p_fue_vec./Cnt;
        P_FUE_Mat_W{i} = 10.^(p_fue_vec./(10*Cnt));
        P_mean_FUE = [P_mean_FUE mean(P_FUE_Mat{i})];
        P_max_FUE = [P_max_FUE max(P_FUE_Mat{i})];
        P_min_FUE = [P_min_FUE min(P_FUE_Mat{i})];
        P_sum_FUE = [P_sum_FUE sum(P_FUE_Mat_W{i})];
    end
end