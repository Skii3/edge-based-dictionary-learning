function [ I2,w ] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin)
% patch_size:提取的patch的大小
% L：patch_group的大小
% step_patch：取patch时的step
% w：存放每个像素的权重矩阵
% L_select_width、L_select_height：在L_select_width * L_select_height范围内选取相似的L个patch
% m:每个patch的长与宽
% I：原始有噪声图片
% t：噪声的估计方差
global D_modify_count 
global D_modify
global D_modify_k

L_select_width_origin = L_select_width;
L_select_height_origin = L_select_height; 
[height, width] = size(I2);
flag = 2;               % 2表示空缺区域补满阶段，3表示空缺区域必须补满阶段
[I2_sort,index] = sort(reshape(I2',1,width * height));                      % 对去噪图像的像素值进行升序排列
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
if I2_sort(ind) ~= 0                                            % 确保ind所指的像素为0
    fprintf('ok');
else                                                                    % 存在像素值为0的点，则需要对其进行单独的svd处理
    
%     while(1)
%         repeat_choice = unidrnd(m);
%         if repeat_choice >= 3
%             break
%         end
%     end
    repeat_choice = 3;
    count_center = 0;
    x = floor(index(ind) / width) + 1;                                       % 计算得到该像素值的坐标位置
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
        %% 得到范围框左上角点的坐标
        x = X(ii);
        y = Y(ii);
        temp_L_height_index = x - round(L_select_height / 3);               % 计算范围框左上角点的纵坐标和横坐标（尽量使点x，y在中间）
        temp_L_width_index = y - round(L_select_width / 3);
        if temp_L_height_index < 1                                          % 若左上角点的纵坐标越界
            temp_L_height_index = 1;
        elseif temp_L_height_index + L_select_height - 1 > height           % 若左上角点的纵坐标造成范围框越界
            temp_L_height_index = height - L_select_height + 1;
        end
        if temp_L_width_index < 1                                          % 若左上角点的横坐标越界
            temp_L_width_index = 1;
        elseif temp_L_width_index + L_select_width - 1 > width           % 若左上角点的横坐标造成范围框越界
            temp_L_width_index = width - L_select_width + 1;
        end
        if L_select_height >= height || L_select_width >= width
            flag = 3;
            L_select_height = height;
            L_select_width = width;
            temp_L_height_index = 1;
            temp_L_width_index = 1;
        end
        %% 得到包含未赋值像素点的中心patch的值
        temp_p_height_index = x - round(m / 3);               % 计算范围框左上角点的纵坐标和横坐标（尽量使点x，y在中间）
        temp_p_width_index = y - round(m / 3);
        if temp_p_height_index < 1                                          % 若左上角点的纵坐标越界
            temp_p_height_index = 1;
        elseif temp_p_height_index + m - 1 > height           % 若左上角点的纵坐标造成范围框越界
            temp_p_height_index = height - m + 1;
        end
        if temp_p_width_index < 1                                          % 若左上角点的横坐标越界
            temp_p_width_index = 1;
        elseif temp_p_width_index + m - 1 > width           % 若左上角点的横坐标造成范围框越界
            temp_p_width_index = width - m + 1;
        end
        center(ii,:) = reshape(I_origin(temp_p_height_index:temp_p_height_index + m - 1,temp_p_width_index:temp_p_width_index + m - 1)',1,m * m);             
    end
    [I2,w] = step1_2_32_2( temp_L_height_index,temp_L_width_index, patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m,I,t,flag,I_origin,center);

    %% 递归调用，直到不存在零像素点
    L_select_width = L_select_width_origin;
    L_select_height = L_select_height_origin;
    [I2,w] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I_origin);

end

end


