function [ I2,w ] = step1_2_32_2( L_height_index,L_width_index, patch_size,L,I2,step_patch,w,L_select_width,L_select_height,m ,I,t,flag,I_origin,center )
% L_height_index：范围框左上角纵向的坐标
% L_width_index：范围框左上角横向的坐标
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
global D_modify_ave
%% step1——patch grouping: block matching method
L_index_begin = [L_height_index L_width_index];                                                     % 取L范围框左上角点的坐标
patch_height_index = [L_index_begin(1):step_patch:L_index_begin(1) + L_select_height - 1 - m, ...    % 范围框内patch左上角初始位置的索引
    L_index_begin(1) + L_select_height - m];
patch_width_index = [L_index_begin(2):step_patch:L_index_begin(2) + L_select_width - 1 - m, ...    % 范围框内patch左上角初始位置的索引
    L_index_begin(2) + L_select_width - m];
ind = 0;
number_patch = length(patch_width_index) * length(patch_height_index);
y = zeros(number_patch, patch_size);
y2 = zeros(number_patch, patch_size);
% P = zeros(patch_size, L + 1);
for k1 = 1:m
    for k2 = 1:m
        ind = ind + 1;
        y(:,ind) = reshape(I(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);                % 得到由范围框内的像素构成的patch
        y2(:,ind) = reshape(I_origin(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);  
    end
end
% 使用欧几里得距离，对中心patch，在范围框内，寻找其L个最近邻
center_patch_num = 1;
repeat_choice = size(center,1);
for iii = 1:center_patch_num^2
    if flag == 2
        for i_center = 1 : repeat_choice
            patch_center = center(i_center,:);
            PC_select2 = normxcorr2(patch_center',D_modify_ave);
            temp_normxcorr(i_center) = max(PC_select2(patch_size,:));
        end
        i_center = find(temp_normxcorr == max(temp_normxcorr));
        patch_center = center(i_center,:);                              % 挑选新的中心区块             
    else
        patch_center = center(floor(repeat_choice / 2) + 1,:);
    end
    
    
    
     % 先归一化之后再计算欧几里得距离，与原文不同
    y2 = y2(:,:) - repmat(mean(y2,2),1,patch_size) * 0.2;
    patch_center = patch_center - repmat(mean(patch_center),1,patch_size) * 0.2;    
    temp_dis = y2(:,:) -  repmat(patch_center,number_patch,1);               % 计算其他的patch到中心patch之间的欧几里得距离                  %%%%
    temp_dis2 = sum(temp_dis.^2,2);
    [temp_dis2_sort,index] = sort(temp_dis2(:));                                         % 按照欧几里得距离进行排序
    L = sum(sum(temp_dis2_sort < t^2 * 100 * 2)) - 1;
    if L <= 1
        L = 2;
    end
    P = zeros(patch_size, L + 1);
    P(:,1:L + 1) = y(index(1:L + 1),:)';                                    % 按距离从近到远取出L个patch
    %% step2——SVD_based denoising
    P_ave_P = mean(P');
    P_SVD = zeros(size(P));
    P_origin = P;
    P = P - repmat(P_ave_P',1,L + 1);
    [U_P,S,V] = svd(P);                                                       % 进行SVD变换
    [M,N] = size(P);
    [n1,n2] = size(S);
    sum_k_1 = 0;
    k = 0;
    for j2 = min(n1,n2) : -1 : 2
        sum_k_1 = sum_k_1 + S(j2,j2)^2;                                     % 对于通过svd变换得到的特征值，从后往前求取平方和
        sum_k = sum_k_1 + S(j2 - 1,j2 - 1)^2;
        if sum_k_1 <= t^2 * (L + 1) * patch_size && sum_k >= t^2 * (L + 1) * patch_size             % 当满足论文所给条件之后
            k = j2 - 1;                                                     % 获得应取的维数
            break;
        end
    end
    
    for i_D = 1 : D_modify_count
        U = reshape(D_modify(i_D,:,1:D_modify_k(i_D)),patch_size,D_modify_k(i_D));               % 得到主成分
        alpha = U' * P;
        P_recon(i_D,:,:) = U * alpha;                % 得到重构矩阵
        recon_error(i_D,1) = sum(sum((reshape(P_recon(i_D,:,:),size(P,1),size(P,2)) - P) .^ 2));
    end
    [recon_error_sort,recon_error_index] = sort(recon_error,'ascend');
    k_select = 0;
    U_index = 1;
    flag_add = 0;
    while(1)
        U_select(:,k_select + 1:k_select + D_modify_k(recon_error_index(U_index))) = ...
            reshape(D_modify(recon_error_index(U_index),:,1:D_modify_k(recon_error_index(U_index))),...
            1,patch_size,D_modify_k(recon_error_index(U_index)));  
        k_select = k_select + D_modify_k(recon_error_index(U_index));
        U_index = U_index + 1;
        q = orth(U_select);
        alpha = q' * P;
        P_recon_qr = q * alpha + repmat(P_ave_P',1,L + 1);                % 得到重构矩阵
        recon_error = sum(sum((P_recon_qr(:,:) - P_origin).^2));
        fprintf('%d,%f\n',U_index,recon_error);
        if recon_error < t^2 * size(P_recon_qr,1) * size(P_recon_qr,2) || U_index == D_modify_count || rank(q) > 10
            if rank(q) > 10
                flag_add = 1;
            end
            break;
        end
    end
    if U_index ~= D_modify_count && flag_add == 0
        k_select = rank(U_select);
        w_temp = 0;
        if k_select < L + 1
            w_temp = 1 - k_select/(L + 1);                                             % 计算该像素的权重
        elseif k_select == L + 1
            w_temp = 1/(L + 1);
        end
        fprintf('%d,%d\n',k_select,L);
%         figure
%         subplot(1,2,1),imdisp(tiledict(P_origin, [10 10]));
%         subplot(1,2,2),imdisp(tiledict(P_recon_qr, [10 10]));
        for ii = 1 : L + 1
            index_height = floor(index(ii) ./ length(patch_width_index)) + 1;
            index_width = index(ii) - length(patch_width_index) * (index_height - 1);
            if index_width == 0
                index_width = length(patch_width_index);
                index_height = index_height - 1;
            end
            X = patch_height_index(index_height);
            Y = patch_width_index(index_width);
            temp = w(X:X + m - 1,Y:Y + m - 1) + w_temp .* ones(m,m);
            if sum(sum(temp == 0)) > 0
                dbstop step1_2_3_2.m at 135
            end     
            I2(X:X + m - 1, Y: Y + m - 1) = ((I2(X:X + m - 1, Y: Y + m - 1) .* w(X:X + m - 1,Y:Y + m - 1) ...
                    + reshape(P_recon_qr(:,ii)',m,m)' * w_temp)) ./ temp;
            w(X:X + m - 1,Y:Y + m - 1) = temp;
        end
        clear P_recon_qr P_origin
    else                                % 添加词典典元
        D_modify_count = D_modify_count + 1;
        D_modify(D_modify_count,:,1:k) = reshape(U_P(:,1:k),1,patch_size,k);
        D_modify_k(D_modify_count,1) = k;
        D_modify_ave(:,D_modify_count) = P_ave_P';
        
        P_SVD(:,:) = U_P(:,1:k) * S(1:k,1:k) * V(:,1:k)' + repmat(P_ave_P',1,L + 1);
        w_temp = 0;
        if k < L + 1
            w_temp = 1 - k/(L + 1);                                             % 计算该像素的权重
        elseif k == L + 1
            w_temp = 1/(L + 1);
        end
%         figure
%         subplot(1,2,1),imdisp(tiledict(P_origin, [10 10]));
%         subplot(1,2,2),imdisp(tiledict(P_recon_qr, [10 10]));
        for ii = 1 : L + 1
            index_height = floor(index(ii) ./ length(patch_width_index)) + 1;
            index_width = index(ii) - length(patch_width_index) * (index_height - 1);
            if index_width == 0
                index_width = length(patch_width_index);
                index_height = index_height - 1;
            end
            X = patch_height_index(index_height);
            Y = patch_width_index(index_width);
            temp = w(X:X + m - 1,Y:Y + m - 1) + w_temp .* ones(m,m);
            if sum(sum(temp == 0)) > 0
                dbstop step1_2_3_2.m at 166
            end         
            I2(X:X + m - 1, Y: Y + m - 1) = ((I2(X:X + m - 1, Y: Y + m - 1) .* w(X:X + m - 1,Y:Y + m - 1) ...
                    + reshape(P_SVD(:,ii)',m,m)' * w_temp)) ./ temp;
            w(X:X + m - 1,Y:Y + m - 1) = temp;
        end  
    end

end
end



