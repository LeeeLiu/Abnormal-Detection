

% ���ߣ�LT
% ���������ã�
	% �춯���⣨LT �� YSH���� �����Ż���

% ���룺
    % SP����������(ԭʼ����)����ȡ SP �����������μ� TXT
	% K1,K2, �춯�����2����ֵ������ �����ʼ����
		% K1,K2 ԽС��Խ���У��������춯��Խ�ࡣ
	% N���춯����ĵ�3����ֵ��
		% N ��ʾ ���� �춯�� �ı�������ʼ��Ϊ1����ʾ ���� ȫ�� �춯�㡣
	
% �����
    % wave_filtered���˲�����(ǰ24��Ϊ0)����������(��һ��)��ʱ���(�ڶ���)��
    % optimal_abnorm���춯�㡾��������(��һ��)��ʱ���(�ڶ���)��	
	% ������ֵ optimal_K1��optimal_K2��optimal_N
	

% ��Ҫ˼��
	% ���裺
		% A Ԥ�����˲�
		% B ��һ�μ�⣨baseline��YSH������
		% C �ڶ��ι��ˣ�LT������
	
	% A Ԥ�����˲�
		% ���� ԭʼ����SP��һ�ײ�֣���������Ԫ�ز�ֵ ����ֵ����
		% ���� ��ֵ ���� ��������*�����ֵ-��Сֵ��
		% ��� SP �У������� С�� �� �� ��ܶ࣬�� �Ӵ��˲�ǿ�ȣ����� ��ֵ�˲� ����

		
	% B ��һ�μ�⣺
	% K1��K2 ��ʼ������С��0.001��0.005����N��ʼ��Ϊ1
		% ���� baseline���� ���еظ���ļ��� ����� �춯�㡣
		% ���� ��Щ �춯�㣬���� �ڶ��ι��ˡ�
    
	% C �ڶ��ι���		
		% ����B���������춯�㣨��ʱK1��K2�ѹ̶��������� �춯�� ������ͬʱ ���� ������״��
		% �����ʧ���ݶ��½� ���� �Ż� N��
			% ���������ж�����˼����YSH�ṩ��
				% �춯�����N , 
				% ���ߣ��˲������ߣ��춯�����ߣ������֮�
			% ���ڴˣ�һ������� 4 ����ʧ������
				% ������ʽ�����ڱ��ļ������� �ؼ��� '% ��ʧ todo'
			
		% 4 ����ʧ ���ܶԱ�
			% ��ʧ��1���ͣ�4��Ŀǰ�ǱȽϺ��õġ�
			% ��ʧ��2���ͣ�3�����ȶ�����������£�Ҫô��һ�ξͽ����ˣ�ҪôN��һֱͣ�ڳ�ʼֵ1��
		
		% �� ����AԤ����֮��� �˲����Σ����� ���ײ�֣��� ��������
			% ȡ�� ���ײ������ ��N %��ȡ ��ģ� ��ʱ��� index��
			% �� �˲������� �ҵ� index ��Ӧ�ĵ㣬
			% ��Щ�㣬���� ��ǰ ���ŵ� �춯�����С�
			% ������������Щ �춯�����У����� ������ optimal_abnorm��

% ����ʼ
    clc;
    %    rng(1);
    rng('shuffle')        % ���ݵ�ǰʱ���ʼ��������
    
    %     ��ȡ MAT
    % load 'D:\2-���빤��\Z-Dataset\01-Internship\BaiDuYun\ÿһ������-MAT\C_latest_parameter.mat'
    % todo 4
    load   'D:\2-���빤��\Z-Dataset\01-Internship\�������ݼ�1\6-parameter.mat';
    total_num = size(SP,1);
    ratio_1 = 0.0;
    ratio_2 = 0.2;
    Start = 1 + ratio_1 *total_num;
    End = round( ratio_2 *total_num );
    len = End - Start +1;
    SP = SP(Start: End, :);
    %     ��ȡ MAT
    figure(1);
    plot(SP) ; % �鿴 �Ƿ� ����
	title('ԭʼ����');

        
    % ���յõ����춯�㣬��ֵ
    % ����ֵ��
    optimal_K1=[];
    optimal_K2=[];
    optimal_N=[];
    optimal_abnorm = [];
    % ����ֵ
    
    % ��ֵ ��ʼ��
    % todo1
    K1 = 0.001 +0.004*rand(1,1);
    K2 = 0.001 +0.004*rand(1,1);
    
    N = 1.0; 
    iteration = 100;
    
    %     K1=0.005; K2=0.05;    
    %     K1 = rand(1,1);  K2 = rand(1,1);    
    % ��ֵ
    
    % Adam����
    beta1 = 0.9;
    beta2 = 0.999;
    eps = 1e-8;
    
    % todo2
    a = 0.01;
%     a = 0.001;
    % Adam����
    
    % ��ʧ
    loss = intmax;
    % ��ʧ
    
    % (filter_num+1)/2 = 25ǰ���Ϊ0
    % �˲�����
    filter_num=49;
    % �˲�����

% 1Ԥ�����˲�
        % 	һ���˲� �õ� SP_1
        SP_1 = movmedian(SP, filter_num);
		% ����matlab�⺯�������� movmedian�����������Ի�ȡ��Ӧ��C����.

	ratio = 0.005;
	delta = ratio *(max(SP_1) - min(SP_1));
	SP_2 = SP_1;
	ratio_2 = 0.001;
	%
    diff_1 = abs( diff(SP_2) );
    index = find(diff_1<delta);     % �ҳ� diff_1 �� С�� delta �� ʱ���index			

	 if  size(index,1)  >= ratio_2 * len
            winSize = round(size(index,1)/80);
			SP_2 = movmedian(SP_2,  winSize );
     end
% 1Ԥ�����˲�

        % 	�õ� �˲��� ���� SP_2
        SP = SP_2;
        figure(2);
        plot(SP);       % �鿴 Ԥ�����˲� Ч��
		title('Ԥ�����˲�');

		clear func_AdamGradientDescent;         %  ���persistent ��̬���������� Ӱ�� ��һ�� ������


        % 2 �����ʼ��K��
        % 		��� �춯������ С�� 0.001*�������� ���� ��ʼ��K��N��
        step = 100;
        for  i = 1 : step
                [~, abnorms] = func_AbnormDetect(SP(:,1),K1,K2,filter_num);    
                if size(abnorms,1) >= 0.001*len
                    break;
                elseif  i==iteration
                    fprintf('��ʼ��ʧ�ܣ�����������');
                    % break;
                    exit(1);
                else
                    K1 = 0.001 +0.099*rand(1,1);
                    K2 = 0.001 +0.099*rand(1,1);
                               N = 1.0;
                end
        end

% % ������
       % ���
       % �˲�+�춯����
       [wave_filtered, abnorms] = func_AbnormDetect(SP(:,1),K1,K2,filter_num);
       wave_filtered = wave_filtered(all(~isnan(wave_filtered),2),:);       % ɾ������NaN����      
       % ������(��һ��)�� ���ײ�֣����� ʱ���(�ڶ���)
       second_order_dif = ( ( abs( diff(  abs(diff( abnorms(:, 1) ) ) ) ) ) )    ;       
       % ƴ�� ʱ��� �����ײ����� ����2����
       second_order_dif = [second_order_dif, abnorms(3:end, 2) ];       
       % �� ���ײ��(��һ��)����С����
       sorted_dif = sortrows( second_order_dif, 1 );

       for  k = 1 : iteration
           % ȡ �� N% �� ��
           cut_point = 1+round( (1-N)*size(sorted_dif, 1) );
           firstN_timestamp = sorted_dif( cut_point: end , 2);
           %  �� ǰһ���� Ҳ �ϱ�
           firstN_timestamp = [firstN_timestamp;firstN_timestamp-1];
           % ȥ��
           firstN_timestamp=unique(firstN_timestamp,'rows');
           % �� wave_filtered �� ȡ�� ʱ��� Ϊ firstN_timestamp �� ��
           firstN_abnorm = wave_filtered( firstN_timestamp+(filter_num+1)/2 , :);
           firstN_abnorm = sortrows(firstN_abnorm, 2);
           
		   % auc �����������
           auc = trapz( wave_filtered(:,2) ,wave_filtered(:,1) );
           auc_hat = trapz( firstN_abnorm(:,2) ,firstN_abnorm(:,1) );
           rectangle =  (eps+(max(wave_filtered(:,1))-min(wave_filtered( (filter_num+1)/2:end ,1) ))*len);
           y        = auc /rectangle;
           y_hat = auc_hat /rectangle;
           
           % ��ʧ todo
           % ��1��MSE  
			% f = @(N)    N + 1e+04*(  y - y_hat  )^2;
           
         % ��2�� ������  �е���  [��һЩ���Σ�������N��һֱ��1.0]
           CE = - y*log(y_hat)  - (1-y)*log(1-y_hat);
           f = @(N)    N * CE  ;
           
           % ��3������log     [��һЩ���Σ�������N��һֱ��1.0]
%               f = @(N) log( eps  + N )/log(0.05) + ( log(  auc / auc_hat  ) )^2  ; 
           % ��4��һ��log
%              f = @(N)  N * ( log(  auc / auc_hat  ) )^2  ;  
     
           % ���
           
           % GD
           fprintf("��%d�Σ�N=%f, K1=%f, K2=%f, f(N)=%f\n", k, N,K1,K2, f(N));
           
           x= N;
           x_grad = (f(x+a) - f(x-a)) / (2*a);
           [x_fin, xn] = func_AdamGradientDescent(f, x,   x_grad,  a,  eps, beta1, beta2);                      
           N  = (xn);
           
           if N<=0 || N>1
               N = 1;
           end                      
           
           % GD
           
           %             if f(N,K1,K2) < loss
           loss = f(N);
           
           optimal_N=N;
           optimal_abnorm =  firstN_abnorm;
          
           my_plot( wave_filtered, optimal_abnorm);
           str = sprintf('N=%f, K1=%f, K2=%f.fig', N,K1,K2);
           path = strcat('�û���-��ʧ1-�м���','\', str);
           savefig('final.fig');
           %                 savefig(path);
           clf(figure(3));
           %             end
            
           if  size( optimal_abnorm,1 )/len   <  1000 / ( 24*60*60 )
               break;
           end           
           
           if x_fin
               fprintf('�ݶ�̫ƽ�����½������ˡ����������������߸������ݡ�');
               break;
           end
           
       end

delete(figure(1));
delete(figure(2));
delete(figure(3));

% % �������
