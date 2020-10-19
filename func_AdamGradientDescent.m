

function [fin, x] = func_AdamGradientDescent(f, x,  grad,  a, eps, beta1, beta2)
% Adam���ݶ��½� ��ֵ�Ż�
% ����	
    % f - һ�����н�����ʽ�ĺ��������� f = @(X) 2*X ;
    % x - ���Ż��ı���
    % a - ����
    % eps - epsilon��һ����С��С��ֵ�����ڷ�ĸ���������0
    % beta1 - Adam����ֵ����
    % beta2 - Adam����ֵ����
% ���
	% x���Ż�һ��֮��ı���
    % fin�����ֵΪ1��˵���ݶ�̫С��ƽ�������޷������Ż�����ֹͣ�����������
    
    persistent k;		% persistent���� �൱�� C/C++ ��̬����
    persistent m;
    persistent v;
    persistent beta1_;
    persistent beta2_;
    
    if isempty(k)
        k = 1;
        m = grad;
        v = grad.^2;
        beta1_ = beta1;
        beta2_ = beta2;
    else
        k = k + 1;
        m = beta1*m + (1-beta1)*grad;
        v = beta2*v + (1-beta2)*grad.^2;
        beta1_ = beta1_ * beta1;
        beta2_ = beta2_ * beta2;
    end
    
    m_hat = m / (1-beta1_);
    v_hat = v / (1-beta2_);
    
    fin = abs(grad) <= 1e-3;
    x = x - a / (sqrt(v_hat) + eps) * m_hat;
end
