

function my_plot( SP_filter, data_sample)
% 输入
	% SP_filter：滤波后 波形
	% data_sample：异动点
	
% 输出	
	% 图figure(3)：异动点 描绘在 滤波后曲线 上
	
    figure(3);
    % 2 滤波后-波形
    plot(SP_filter(:,end), SP_filter(:,1));
    hold on
    % 3 异动点
    plot(data_sample(:,end), data_sample(:,1), 'o' , 'Color',[1 0 0]);
    legend( '滤波后波形', '异动点');

    xlabel('index(时间)');
    ylabel('特征值(A相电流)');
    title('异动点检测');

end
