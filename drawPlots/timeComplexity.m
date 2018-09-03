
clear;
clc;
%%
T0 = 'Aug26/T1/pro_x2_IL_%d_%d.mat';
T1 = 'Aug26/T1/pro_x3_IL_%d_%d.mat';
T2 = 'Sep2/T1/pro_x2_CL_%d_%d.mat';
T3 = 'Sep2/T1/pro_x3_CL_%d_%d.mat';
% T4 = 'Aug23/T2/pro_greedy_%d_%d.mat';
%%
T_vec0 = TimeComplexity(T0);
T_vec1 = TimeComplexity(T1);
T_vec2 = TimeComplexity(T2);
T_vec3 = TimeComplexity(T3);
% T_vec4 = TimeComplexity(T4);
%%
figure;
hold on;
grid on;
box on;
plot(T_vec0, '--or', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r');%, 'MarkerEdgeColor','b');
plot(T_vec1, '--ob', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b');%, 'MarkerEdgeColor','b');
plot(T_vec2, '--dr', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','r');%, 'MarkerEdgeColor','b');
plot(T_vec3, '--db', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','b');%, 'MarkerEdgeColor','b');
% plot(T_vec4, '--*k', 'LineWidth',1.3,'MarkerSize',8, 'MarkerFaceColor','k');
% title('Time duration of proposed RF','FontSize',14, 'FontWeight','bold');
xlabel('FBS Numbers','FontSize',12);
ylabel('seconds','FontSize',12);
xlim([1 10]);
ylim([10 20]);
legend({'IL+$\mathcal{X}_1$','IL+$\mathcal{X}_2$', 'CL+$\mathcal{X}_1$', 'CL+$\mathcal{X}_2$'},'Interpreter','latex','FontSize',12);
% legend({'proposed RF','[9]'},'FontSize',14, 'FontWeight','bold');

%%
function T_vec = TimeComplexity(input)
    T_vec = [];
    for i=1:10
        fprintf('FBS num = %d\t', i);
        Cnt = 0;
        T = 0.0;

        for j=1:500
            s = sprintf(input,i,j);
            filename = strcat(s);
            if exist(s)
                load(filename);
                T = T + QFinal.time;
                Cnt = Cnt+1;
            end
        end
        fprintf('Total Cnt = %d\n',Cnt);
        T_vec = [T_vec T/(i*Cnt)];
    end
end
