function [ I2,w ] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin)
% patch_size:��ȡ��patch�Ĵ�С
% L��patch_group�Ĵ�С
% step_patch��ȡpatchʱ��step
% w�����ÿ�����ص�Ȩ�ؾ���
% L_select_width��L_select_height����L_select_width * L_select_height��Χ��ѡȡ���Ƶ�L��patch
% m:ÿ��patch�ĳ����
% I��ԭʼ������ͼƬ
% t�������Ĺ��Ʒ���
global D_modify_count 
global D_modify
global D_modify_k

L_select_width_origin = L_select_width;
L_select_height_origin = L_select_height; 
[height, width] = size(I2);
flag = 2;               % 2��ʾ��ȱ�������׶Σ�3��ʾ��ȱ������벹���׶�
[I2_sort,index] = sort(reshape(I2',1,width * height));                      % ��ȥ��ͼ�������ֵ������������
ind = 1;
while(1)
    if I2_sort(ind) < 0
        ind = ind + 1;
    else
        break;
    end
end
ind_0_begin = ind;
ind = 1;
while(1)
    if I2_sort(ind) <= 0
        ind = ind + 1;
    else
        break;
    end
end
ind_0_end = ind - 1;
fprintf('ind_0_end:%d\n',ind_0_end);
% for iii = ind_0_begin:1:ind_0_end
ind = ind_0_begin;
if I2_sort(ind) ~= 0                                            % ȷ��ind��ָ������Ϊ0
    fprintf('ok');
else                                                                    % ��������ֵΪ0�ĵ㣬����Ҫ������е�����svd����
    
%     while(1)
%         repeat_choice = unidrnd(m);
%         if repeat_choice >= 3
%             break
%         end
%     end
    repeat_choice = 3;
    count_center = 0;
    x = floor(index(ind) / width) + 1;                                       % ����õ�������ֵ������λ��
    y = index(ind) - width * (x - 1);
    if y == 0        
        x = x - 1;
        y = width;
    end
    X = [x - (repeat_choice - 1) : 1 : x + (repeat_choice - 1)];
    Y = [y - (repeat_choice - 1) : 1 : y + (repeat_choice - 1)];
    if X(1) < 1
        flag = 3;
        X = X + ones(size(X)) * (abs(X(1)) + 1);
    end
    if Y(1) < 1
        flag = 3;
        Y = Y + ones(size(Y)) * (abs(Y(1)) + 1);
    end
    if X((repeat_choice - 1) * 2 + 1) > size(I,1)
        flag = 3;
        X = X - ones(size(X)) * (abs(X((repeat_choice - 1) * 2 + 1) - size(I,1)));
    end
    if Y((repeat_choice - 1) * 2 + 1) > size(I,1)
        flag = 3;
        Y = Y - ones(size(Y)) * (abs(Y((repeat_choice - 1) * 2 + 1) - size(I,2)));
    end
    
    for ii = 1 : length(X)
        %% �õ���Χ�����Ͻǵ������
        x = X(ii);
        y = Y(ii);
        temp_L_height_index = x - round(L_select_height / 3);               % ���㷶Χ�����Ͻǵ��������ͺ����꣨����ʹ��x��y���м䣩
        temp_L_width_index = y - round(L_select_width / 3);
        if temp_L_height_index < 1                                          % �����Ͻǵ��������Խ��
            temp_L_height_index = 1;
        elseif temp_L_height_index + L_select_height - 1 > height           % �����Ͻǵ����������ɷ�Χ��Խ��
            temp_L_height_index = height - L_select_height + 1;
        end
        if temp_L_width_index < 1                                          % �����Ͻǵ�ĺ�����Խ��
            temp_L_width_index = 1;
        elseif temp_L_width_index + L_select_width - 1 > width           % �����Ͻǵ�ĺ�������ɷ�Χ��Խ��
            temp_L_width_index = width - L_select_width + 1;
        end
        if L_select_height >= height || L_select_width >= width
            flag = 3;
            L_select_height = height;
            L_select_width = width;
            temp_L_height_index = 1;
            temp_L_width_index = 1;
        end
        %% �õ�����δ��ֵ���ص������patch��ֵ
        temp_p_height_index = x - round(m / 3);               % ���㷶Χ�����Ͻǵ��������ͺ����꣨����ʹ��x��y���м䣩
        temp_p_width_index = y - round(m / 3);
        if temp_p_height_index < 1                                          % �����Ͻǵ��������Խ��
            temp_p_height_index = 1;
        elseif temp_p_height_index + m - 1 > height           % �����Ͻǵ����������ɷ�Χ��Խ��
            temp_p_height_index = height - m + 1;
        end
        if temp_p_width_index < 1                                          % �����Ͻǵ�ĺ�����Խ��
            temp_p_width_index = 1;
        elseif temp_p_width_index + m - 1 > width           % �����Ͻǵ�ĺ�������ɷ�Χ��Խ��
            temp_p_width_index = width - m + 1;
        end
        center(ii,:) = reshape(I_origin(temp_p_height_index:temp_p_height_index + m - 1,temp_p_width_index:temp_p_width_index + m - 1)',1,m * m);             
    end
    [I2,w] = step1_2_32_2( temp_L_height_index,temp_L_width_index, patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m,I,t,flag,I_origin,center);

    %% �ݹ���ã�ֱ�������������ص�
    L_select_width = L_select_width_origin;
    L_select_height = L_select_height_origin;
    [I2,w] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin);

end

end


