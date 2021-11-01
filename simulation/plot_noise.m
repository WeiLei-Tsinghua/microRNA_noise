clc
clear

%% set the kind of noise
% 0: all noise; i: noise contributed by reaction i
reaction_i = 0;

%% load data and extract related parameters
load exp_CV.mat
load exp_mean.mat

data_size = size(exp_CV);
N_alter = data_size(1);

%% plot
exp_CV_plot = exp_CV(:, :, reaction_i + 1);

subplot(1,2,1)

for i = 1:N_alter
    Ctrl = exp_mean(1, :);
    FC = log10(exp_mean(i,:)./Ctrl);
    semilogx(exp_mean(i,:), FC, 'LineWidth', 1);
    legend_str{i} = ['condition' num2str(i)];
    hold on
end

title('Fold change')
xlabel('log10(expression) (a.u.)')
ylabel('log10(fold change)')
legend(legend_str)

%% noise
subplot(1,2,2)
for i = 1:N_alter
    semilogx(exp_mean(i,:), exp_CV_plot(i,:), 'LineWidth', 1);
    legend_str{i} = ['condition' num2str(i)];
    hold on
end

title('CV')
xlabel('log10(expression) (a.u.)')
ylabel('CV')
legend(legend_str)
