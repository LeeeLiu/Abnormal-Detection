

% （对比 Baseline）
%   调用 YSH的方法
[wave_filtered, abnorms] = func_AbnormDetect(SP,0.005,0.05,filter_num);
my_plot( wave_filtered, abnorms);
savefig('对比-数据XXX-K1=0.005, K2=0.05.fig');

