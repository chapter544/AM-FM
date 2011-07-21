function imstat(x)

min_x = min(x(:));
max_x = max(x(:));
var_x = var(x(:));
std_x = std(x(:));
mean_x = mean(x(:));

display(sprintf('Min : %f', min_x));
display(sprintf('Max: %f', max_x));
display(sprintf('Mean: %f', mean_x));
display(sprintf('Variance : %f', var_x));
display(sprintf('Std: %f', std_x));
