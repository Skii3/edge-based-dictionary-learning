% 论文：An efficient SVD_based method for image denoising 基于边缘的加速
close all
clear
clear global
%% read the noise image
global D_modify_count 
global D_modify
global D_modify_ave
global D_modify_k
D_modify_count = 0;
% I = imread('lena.jpg');                                           % 读取图片
I = imread('barbara.png'); 
I_origin = double(I);
I_size = size(I);                                                     % 得到通道数，若为3通道则转化为1通道，若为1通道则不处理
if length(I_size) == 3 && I_size(3) == 3
    I = rgb2gray(I);
end
imwrite(I,'原图像.jpg');                                            % 保存原始图像
imshow(I),title('原始图像');                                        % 显示原始图像
t = 30;
randn('seed', 0);                                                   % 生成随机种子
noise = t * randn(size(I));                                        % 通过随机种子得到噪声图像
I = double(I);
I = I + noise;                                                      % 得到加入噪声之后的图
I = imresize(I, 1);
figure,imshow(I / 255),title('加高斯白噪声');                        % 保存加入噪声之后的图
imwrite(I / 255,'加高斯白噪声.jpg');
tic
%% estimate noise standard devition by computing median absolute devition of the finest wavelet coefficients

[~,B,C,D] = dwt2(I,'db1');                                          % 二维小波变换
t = median(median(median(abs([B,C,D])))) / 0.6745 ;                  % 估计噪声的标准偏差

f = fspecial('gaussian',[2 2],5);    
I_origin_est = imfilter(I,f,'same');

figure,imshow(I_origin_est,[]);
imwrite(I_origin_est / 255,'g.jpg');

H = fspecial('log',[3 3],10);
edge = imfilter(I_origin_est,H,'replicate');
figure,imshow(edge,[]);
imwrite(edge * 2/(max(max(edge)) - min(min(edge))),'e.jpg');

f2 = fspecial('gaussian',[2 2],5);    
N_origin_est = imfilter(noise,f2,'same');

H2 = fspecial('log',[3 3],10);
edge2 = imfilter(N_origin_est,H,'replicate');
figure,imshow(edge2,[]);

                                        % 二维小波变换
t2 = var(edge2(:));                % 估计噪声的标准偏差

tic
% L_select_width = round(temp(1) / 13.5);                                      % 在L_select_width*L_select_height范围内选取相似的L个patch
% L_select_height = round(temp(2) / 13.5);
% %% 得到亮度最大的block，即边缘最多的部分
% step_L_select_width = 5;                                                    % 取L时范围框的横向移动步长
% step_L_select_height = 5;                                                   % 取L时范围框的纵向移动步长
% L_width_index = [1:step_L_select_width:temp(1) - L_select_width, temp(1) - L_select_width + 1];             % 计算范围框左上角横向的坐标
% L_height_index = [1:step_L_select_height:temp(2) - L_select_height, temp(2) - L_select_height + 1];       % 计算范围框左上角纵向的坐标
% ind = 0;
% number_block = length(L_width_index) * length(L_height_index);
% for k1 = 1:L_select_width
%     for k2 = 1:L_select_height
%         ind = ind + 1;
%         y(:,ind) = reshape(edge(L_width_index + k1 - 1,L_height_index + k2 - 1)',1,number_block);                % 得到由block内的像素构成的patch
%     end
% end
% block_sum = var(y');
% [block_sum_sort,block_sum_index] = sort(block_sum,'descend');
% 
% m = 10;
% patch_size = m * m;
% 
% count = 0;
% for i = 1:number_block
%     index = block_sum_index(i);
%     center_block_x = floor(index ./ length(L_width_index)) + 1;
%     center_block_y = index - length(L_width_index) * (center_block_x - 1);
%     if center_block_y == 0
%         center_block_y = length(L_width_index);
%         center_block_x = center_block_x - 1;
%     end
%     X = L_height_index(center_block_x);
%     Y = L_width_index(center_block_y);
%     flag = 0;
%     for j = 1:count
%         if abs(X - XY_origin(j,1)) < L_select_width / 2 || abs(Y - XY_origin(j,2)) < L_select_height / 2
%             flag = 1;
%             break;
%         end
%     end
%     if flag == 1
%         continue;
%     end
%     
%     count = count + 1;
%     XY_origin(count,1) = X;
%     XY_origin(count,2) = Y;
%     hFig = figure;
%     hAx  = axes;
%     imshow(edge,[],'Parent', hAx);
%     imrect(hAx, [Y, X, L_select_width, L_select_width]);
    w = zeros(size(I));
    I2 = zeros(size(I));
    thr_norm = 0.65;                     % 根据相关系数选择patch group的阈值
    thr_edge_select = t2 * 4;           % 根据edge图中patch的方差挑选中心patch
    %% 在每个block中寻找边缘最多的patch
    step_patch = 1;
    m = 10;
    patch_size = m * m;
    L_index_begin = [1 1];                                                     % 取L范围框左上角点的坐标
    L_select_height = I_size(1);
    L_select_width = I_size(2);
    patch_height_index = [L_index_begin(1):step_patch:L_index_begin(1) + L_select_height - 1 - m, ...    % 范围框内patch左上角初始位置的索引
        L_index_begin(1) + L_select_height - m];
    patch_width_index = [L_index_begin(2):step_patch:L_index_begin(2) + L_select_width - 1 - m, ...    % 范围框内patch左上角初始位置的索引
        L_index_begin(2) + L_select_width - m];

    number_patch = length(patch_width_index) * length(patch_height_index);
    y = zeros(number_patch, patch_size);
    y2 = zeros(number_patch, patch_size);
    
    ind = 0;
    for k1 = 1:m
        for k2 = 1:m
            ind = ind + 1;
            y(:,ind) = reshape(edge(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);                % 得到由block内的像素构成的patch
            y2(:,ind) = reshape(I(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);  
        end
    end
    patch_sum = var(y');
    [patch_sum_sort,patch_sum_index] = sort(patch_sum,'descend');
    % 将边缘最密集的patch作为中心patch
%     ii = patch_sum_index(1);
    count = 0;
    for i = 1:number_patch
        index = patch_sum_index(i);
        if patch_sum(index) < thr_edge_select              % 若边缘密集区域方差小于噪声方差
            break
        end
        center_patch_x = floor(index ./ length(patch_width_index)) + 1;
        center_patch_y = index - length(patch_width_index) * (center_patch_x - 1);
        if center_patch_y == 0
            center_patch_y = length(patch_width_index);
            center_patch_x = center_patch_x - 1;
        end
        X = patch_height_index(center_patch_x);
        Y = patch_width_index(center_patch_y);

        flag = 0;
        for j = 1:count
            if abs(X - XY_origin(j,1)) < m / 2 && abs(Y - XY_origin(j,2)) < m / 2
                flag = 1;
                break;
            end
            center_patch_temp_w = w(X:X + m -1,Y:Y + m - 1);
            center_patch_temp_w_sum = sum(sum(center_patch_temp_w));
            if center_patch_temp_w_sum > 2 * patch_size
                flag = 1;
                break;
            end
        end
        if flag == 1
            continue;
        end

        count = count + 1;
        XY_origin(count,1) = X;
        XY_origin(count,2) = Y;                 % 获得中心patch的位置坐标
        
%         hFig = figure;
%         hAx  = axes;
%         imshow(edge,[],'Parent', hAx);
%         imrect(hAx, [Y, X, m, m]);
        
%         center_patch = y(index,:);
%         PC_select = normxcorr2(center_patch,y);
        %  figure, surf(PC_select), shading flat
        
        center_patch2 = I(X:X + m -1,Y:Y + m - 1);
        PC_select2 = normxcorr2(center_patch2,I);
        
        [P_vec_sort2,P_vec_index2] = sort(PC_select2(:),'descend');
        P_select_num = sum(P_vec_sort2 > thr_norm);
        j_count = 0;
        for j = 1:P_select_num
            [ypeak, xpeak] = find(PC_select2 == P_vec_sort2(j));
            yy = ypeak - m + 1;
            xx = xpeak - m + 1;
            if yy <= 0 || xx <= 0 || yy + m - 1 > I_size(1) || xx + m - 1 > I_size(2)
                continue
            end
            j_count = j_count + 1;
            yoffSet(j_count) = yy;
            xoffSet(j_count) = xx;
%             fprintf('%f,%f\n',yoffSet(j_count),xoffSet(j_count))
            P_group(j_count,:) = reshape(I(yy:yy + m - 1, xx:xx + m - 1),1,m * m);
%             hFig2 = figure;
%             hAx2  = axes;
%             imshow(	I,[],'Parent', hAx2);
%             imrect(hAx2, [xoffSet(j), yoffSet(j), m, m]);
        end
        
%         P_vec = PC_select(:,patch_size);
%         [P_vec_sort,P_vec_index] = sort(P_vec(:),'descend');
%         
%         P_select_num = sum(P_vec_sort > 0.8);
%         P_group = y2(P_vec_sort > 0.8,:);                % 得到patchgroup
%         index2 = index(P_vec_sort > 0.8);               % 得到每个patch对应的索引
        if j_count == 1
            continue;
        end
        P_ave = mean(P_group)';
        L = j_count - 1;
        P = P_group' - repmat(P_ave,1,L + 1); 
        [U,S,V] = svd(P);                                                       % 进行SVD变换
        [M,N] = size(P);
        [n1,n2] = size(S);
        sum_k_1 = 0;
        k = 0;
        for j2 = min(n1,n2) : -1 : 2
            sum_k_1 = sum_k_1 + S(j2,j2)^2;                                                             % 对于通过svd变换得到的特征值，从后往前求取平方和
            sum_k = sum_k_1 + S(j2 - 1,j2 - 1)^2;
            if sum_k_1 <= t^2 * (L + 1) * patch_size && sum_k >= t^2 * (L + 1) * patch_size             % 当满足论文所给条件之后
                k = j2 - 1;                                                                             % 获得应取的维数
                break;
            end
        end
        fprintf('%d,%d\n',k,L);
        P_SVD(:,:) = U(:,1:k) * S(1:k,1:k) * V(:,1:k)' + repmat(P_ave,1,L + 1);
%         if k ~= 0
            D_modify_count = D_modify_count + 1;
            D_modify(D_modify_count,:,1:k) = reshape(U(:,1:k),1,patch_size,k);
            D_modify_k(D_modify_count,1) = k;
            D_modify_ave(:,D_modify_count) = P_ave;
%         end
        % aggregate pixel
        w_temp = 0;
        if k < L + 1
            w_temp = 1 - k/(L + 1);                                             % 计算该像素的权重
        elseif k == L + 1
            w_temp = 1/(L + 1);
        end
        for ii = 1 : L + 1
            X = yoffSet(ii);
            Y = xoffSet(ii);
            temp = w(X:X + m - 1,Y:Y + m - 1) + w_temp .* ones(m,m);
            I2(X:X + m - 1, Y: Y + m - 1) = ((I2(X:X + m - 1, Y: Y + m - 1) .* w(X:X + m - 1,Y:Y + m - 1) ...
                    + reshape(P_SVD(:,ii)',m,m) * w_temp)) ./ temp;
            w(X:X + m - 1,Y:Y + m - 1) = temp;
        end    
        clear P_group yoffSet xoffSet P_SVD
    end
    figure,imshow(I2,[])
    figure,imshow(w,[])
    
    step_patch = 1;
    L = 85;
    L_select_width = round(sqrt(L) * 4.1);                                      % 在L_select_width*L_select_height范围内选取相似的L个patch
    L_select_height = round(sqrt(L) * 4.1);

    [I2,w] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I);
toc
    figure,imshow(I2,[])
    figure,imshow(w,[])

    K = I2;
    MSE = sum(sum((I_origin - K).^2)) / I_size(1) / I_size(2);
    PSNR2 = 20 * log10(255 / sqrt(MSE))



save D_modify2.mat D_modify D_modify_ave D_modify_count D_modify_k