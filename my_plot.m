

function my_plot( SP_filter, data_sample)
% ����
	% SP_filter���˲��� ����
	% data_sample���춯��
	
% ���	
	% ͼfigure(3)���춯�� ����� �˲������� ��
	
    figure(3);
    % 2 �˲���-����
    plot(SP_filter(:,end), SP_filter(:,1));
    hold on
    % 3 �춯��
    plot(data_sample(:,end), data_sample(:,1), 'o' , 'Color',[1 0 0]);
    legend( '�˲�����', '�춯��');

    xlabel('index(ʱ��)');
    ylabel('����ֵ(A�����)');
    title('�춯����');

end
