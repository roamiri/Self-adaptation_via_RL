
clear;
clc;
%%
s = sprintf('oct10/FUE_R_18_time.mat');
filename = strcat(s);
load(filename);

s = sprintf('oct10/FUE_Rref_noshare.mat');
filename = strcat(s);
load(filename);

%%
% figure;
hold on;
grid on;
box on;

fariness = zeros(1,10);
fariness_ref = zeros(1,10);
for i=1:10
vec = C_FUE_Mat_2{i};
% vec_ref = C_FUE_Mat_ref{i};
num = 0.0;
num_ref = 0.0;
denom = 0.0;
denom_ref = 0.0;
n = size(vec,2);
for j=1:n
    num = num + vec(j);
%     num_ref = num_ref + vec_ref(j);
    denom = denom + vec(j)^2;
%     denom_ref = denom_ref + vec_ref(j)^2;
end
    fariness(i) = (num^2)/(n*denom);
%     fairness_ref(i) = (num_ref^2)/(n*denom_ref);
end
plot(fariness, '--og', 'LineWidth',1.2,'MarkerSize',8, 'MarkerFaceColor','g', 'MarkerEdgeColor','b');
% plot(fairness_ref, 'b--.', 'LineWidth',1,'MarkerSize',10);
xlim([1 10]);
ylim([0 1.05]);
% title('Fairness index','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);%, 'FontWeight','bold');
ylabel('Jains Index For Fairness','FontSize',12);%, 'FontWeight','bold');
legend({'IL','CL', '$\rho$'},'Interpreter','latex','FontSize',12);
% legend({'X_1','X_2', 'X_3', 'X_4'},'FontSize',14, 'FontWeight','bold');
% legend({'IL+$\mathbf{X}_2$','CL+$\mathbf{X}_3$'},'FontSize',14, 'FontWeight','bold','Interpreter','latex');