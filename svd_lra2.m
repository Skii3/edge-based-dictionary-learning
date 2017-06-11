% ���ģ�An efficient SVD_based method for image denoising ���ڱ�Ե�ļ���
close all
clear
clear global
%% read the noise image
global D_modify_count 
global D_modify
global D_modify_ave
global D_modify_k
D_modify_count = 0;
% I = imread('lena.jpg');                                           % ��ȡͼƬ
I = imread('barbara.png'); 
I_origin = double(I);
I_size = size(I);                                                     % �õ�ͨ��������Ϊ3ͨ����ת��Ϊ1ͨ������Ϊ1ͨ���򲻴���
if length(I_size) == 3 && I_size(3) == 3
    I = rgb2gray(I);
end
imwrite(I,'ԭͼ��.jpg');                                            % ����ԭʼͼ��
imshow(I),title('ԭʼͼ��');                                        % ��ʾԭʼͼ��
t = 30;
randn('seed', 0);                                                   % �����������
noise = t * randn(size(I));                                        % ͨ��������ӵõ�����ͼ��
I = double(I);
I = I + noise;                                                      % �õ���������֮���ͼ
I = imresize(I, 1);
figure,imshow(I / 255),title('�Ӹ�˹������');                        % �����������֮���ͼ
imwrite(I / 255,'�Ӹ�˹������.jpg');
tic
%% estimate noise standard devition by computing median absolute devition of the finest wavelet coefficients

[~,B,C,D] = dwt2(I,'db1');                                          % ��άС���任
t = median(median(median(abs([B,C,D])))) / 0.6745 ;                  % ���������ı�׼ƫ��

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

                                        % ��άС���任
t2 = var(edge2(:));                % ���������ı�׼ƫ��

tic
% L_select_width = round(temp(1) / 13.5);                                      % ��L_select_width*L_select_height��Χ��ѡȡ���Ƶ�L��patch
% L_select_height = round(temp(2) / 13.5);
% %% �õ���������block������Ե���Ĳ���
% step_L_select_width = 5;                                                    % ȡLʱ��Χ��ĺ����ƶ�����
% step_L_select_height = 5;                                                   % ȡLʱ��Χ��������ƶ�����
% L_width_index = [1:step_L_select_width:temp(1) - L_select_width, temp(1) - L_select_width + 1];             % ���㷶Χ�����ϽǺ��������
% L_height_index = [1:step_L_select_height:temp(2) - L_select_height, temp(2) - L_select_height + 1];       % ���㷶Χ�����Ͻ����������
% ind = 0;
% number_block = length(L_width_index) * length(L_height_index);
% for k1 = 1:L_select_width
%     for k2 = 1:L_select_height
%         ind = ind + 1;
%         y(:,ind) = reshape(edge(L_width_index + k1 - 1,L_height_index + k2 - 1)',1,number_block);                % �õ���block�ڵ����ع��ɵ�patch
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
    thr_norm = 0.65;                     % �������ϵ��ѡ��patch group����ֵ
    thr_edge_select = t2 * 4;           % ����edgeͼ��patch�ķ�����ѡ����patch
    %% ��ÿ��block��Ѱ�ұ�Ե����patch
    step_patch = 1;
    m = 10;
    patch_size = m * m;
    L_index_begin = [1 1];                                                     % ȡL��Χ�����Ͻǵ������
    L_select_height = I_size(1);
    L_select_width = I_size(2);
    patch_height_index = [L_index_begin(1):step_patch:L_index_begin(1) + L_select_height - 1 - m, ...    % ��Χ����patch���Ͻǳ�ʼλ�õ�����
        L_index_begin(1) + L_select_height - m];
    patch_width_index = [L_index_begin(2):step_patch:L_index_begin(2) + L_select_width - 1 - m, ...    % ��Χ����patch���Ͻǳ�ʼλ�õ�����
        L_index_begin(2) + L_select_width - m];

    number_patch = length(patch_width_index) * length(patch_height_index);
    y = zeros(number_patch, patch_size);
    y2 = zeros(number_patch, patch_size);
    
    ind = 0;
    for k1 = 1:m
        for k2 = 1:m
            ind = ind + 1;
            y(:,ind) = reshape(edge(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);                % �õ���block�ڵ����ع��ɵ�patch
            y2(:,ind) = reshape(I(patch_height_index + k1 - 1,patch_width_index + k2 - 1)',1,number_patch);  
        end
    end
    patch_sum = var(y');
    [patch_sum_sort,patch_sum_index] = sort(patch_sum,'descend');
    % ����Ե���ܼ���patch��Ϊ����patch
%     ii = patch_sum_index(1);
    count = 0;
    for i = 1:number_patch
        index = patch_sum_index(i);
        if patch_sum(index) < thr_edge_select              % ����Ե�ܼ����򷽲�С����������
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
        XY_origin(count,2) = Y;                 % �������patch��λ������
        
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
%         P_group = y2(P_vec_sort > 0.8,:);                % �õ�patchgroup
%         index2 = index(P_vec_sort > 0.8);               % �õ�ÿ��patch��Ӧ������
        if j_count == 1
            continue;
        end
        P_ave = mean(P_group)';
        L = j_count - 1;
        P = P_group' - repmat(P_ave,1,L + 1); 
        [U,S,V] = svd(P);                                                       % ����SVD�任
        [M,N] = size(P);
        [n1,n2] = size(S);
        sum_k_1 = 0;
        k = 0;
        for j2 = min(n1,n2) : -1 : 2
            sum_k_1 = sum_k_1 + S(j2,j2)^2;                                                             % ����ͨ��svd�任�õ�������ֵ���Ӻ���ǰ��ȡƽ����
            sum_k = sum_k_1 + S(j2 - 1,j2 - 1)^2;
            if sum_k_1 <= t^2 * (L + 1) * patch_size && sum_k >= t^2 * (L + 1) * patch_size             % ������������������֮��
                k = j2 - 1;                                                                             % ���Ӧȡ��ά��
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
            w_temp = 1 - k/(L + 1);                                             % ��������ص�Ȩ��
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
    L_select_width = round(sqrt(L) * 4.1);                                      % ��L_select_width*L_select_height��Χ��ѡȡ���Ƶ�L��patch
    L_select_height = round(sqrt(L) * 4.1);

    [I2,w] = svd_again2( I2,patch_size,L, step_patch,w,L_select_width,L_select_height,m,I,t,I);
toc
    figure,imshow(I2,[])
    figure,imshow(w,[])

    K = I2;
    MSE = sum(sum((I_origin - K).^2)) / I_size(1) / I_size(2);
    PSNR2 = 20 * log10(255 / sqrt(MSE))



save D_modify2.mat D_modify D_modify_ave D_modify_count D_modify_k