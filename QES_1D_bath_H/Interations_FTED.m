function Ih2=Interations_FTED(L)
% Ih2(��һ����������λ�ã��ڶ�����������λ��)
% Ҫ���һ���������С�ڵڶ�������

% �������
% Ih2.Hp=... 
%        [2,6;
%         5,6;
%         6,7;
%         5,12;
%         7,10;
%         11,12;
%         10,11;
%         11,15];
Ih2.Hp=[(2:L-2).',(3:L-1).'];

% ����bath���
% ��һ������ָ���ʾ����sites��λ�ñ��
% ������ָ���ʾbath��ʹ�õڼ���v
% ���ĸ�ָ��Ϊ0��������ָ����ǰ��Ϊ1��������ָ���ں�
% Ih2.Hb=...        
%       [1,2,1,1;
%        2,3,2,0;
%        4,5,1,1;
%        7,8,2,0;
%        12,13,4,0;
%        9,10,3,1;
%        14,15,4,1;
%        15,16,3,0];
Ih2.Hb=[1,2,1,1;
                L-1,L,2,0];
end