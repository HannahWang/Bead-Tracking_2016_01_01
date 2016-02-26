function strainEnergyStatistics( seres )

filter_se_min = 0;
filter_se_max = 0.000025;
%filter_se_min = 0.000003;
%filter_se_max = 0.00001;
n_plot_per_row = 6;
n_hist_bin = 50;

%
% SE density projection figure
%
Ugs=smooth3(seres.Ucdg_near,'box',[5 5 5]);
Ugs_flatten = reshape(Ugs, prod(size(Ugs)), 1);
max_se = min(max(Ugs_flatten), filter_se_max);
min_se = max(min(Ugs_flatten), filter_se_min);
step = (max_se-min_se)/n_hist_bin;
se_slice = min_se:step:max_se;

% all
figure();
[nelements, centers] = hist(Ugs_flatten(Ugs_flatten>=filter_se_min&Ugs_flatten<=filter_se_max), se_slice);
hist(Ugs_flatten(Ugs_flatten>=filter_se_min&Ugs_flatten<=filter_se_max), se_slice);
xlim([min_se, max_se]);
sprintf('max=%f, median=%f, min=%f', max(Ugs_flatten), median(Ugs_flatten), min(Ugs_flatten))
save('histdata.txt', 'nelements', 'centers', '-ascii')


% sliced histogram for X
figure();
mat_size = size(Ugs);
n_slice = mat_size(2)
for x=1:n_slice
    subplot(ceil(n_slice/n_plot_per_row), n_plot_per_row, x);
    Ugs_filtered = Ugs(:, x, :);
    Ugs_flatten = reshape(Ugs_filtered, prod(size(Ugs_filtered)), 1);
    Ugs_flatten = Ugs_flatten(Ugs_flatten>=filter_se_min&Ugs_flatten<=filter_se_max);
    hist(Ugs_flatten, se_slice);
    xlim([min_se, max_se]);
    title(sprintf('x=%f', seres.xg(1, x, 1)));
    %sprintf('max=%f, median=%f, min=%f', max(Ugs_flatten), median(Ugs_flatten), min(Ugs_flatten))    
end

% color map
for x=1:n_slice
    size(Ugs(:,x,:))
    hh = reshape(Ugs(:,x,:), 47,60);
    HeatMap(hh);
end

% sliced histogram for Y
figure();
mat_size = size(Ugs);
n_slice = mat_size(1)
for y=1:n_slice
    subplot(ceil(n_slice/n_plot_per_row), n_plot_per_row, y);
    Ugs_filtered = Ugs(y, :, :);
    Ugs_flatten = reshape(Ugs_filtered, prod(size(Ugs_filtered)), 1);
    Ugs_flatten = Ugs_flatten(Ugs_flatten>=filter_se_min&Ugs_flatten<=filter_se_max);    
    hist(Ugs_flatten, se_slice);
    xlim([min_se, max_se]);
    title(sprintf('y=%f', seres.yg(y, 1, 1)));    
    %sprintf('max=%f, median=%f, min=%f', max(Ugs_flatten), median(Ugs_flatten), min(Ugs_flatten))
    
end

% sliced histogram for Z
figure();
mat_size = size(Ugs);
n_slice = mat_size(3)
for z=1:n_slice
    subplot(ceil(n_slice/n_plot_per_row), n_plot_per_row, z);
    Ugs_filtered = Ugs(:, :, z);
    Ugs_flatten = reshape(Ugs_filtered, prod(size(Ugs_filtered)), 1);
    Ugs_flatten = Ugs_flatten(Ugs_flatten>=filter_se_min&Ugs_flatten<=filter_se_max);    
    hist(Ugs_flatten, se_slice);
    xlim([min_se, max_se]);
    title(sprintf('z=%f', seres.zg(1, 1, z)));            
    %sprintf('max=%f, median=%f, min=%f', max(Ugs_flatten), median(Ugs_flatten), min(Ugs_flatten))
    
end


end

