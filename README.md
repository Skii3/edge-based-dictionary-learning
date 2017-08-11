# edge-based-dictionary-learning
edge based dictionary learning for image denoising

1. 对有噪图像使用LoG边缘检测检测出边缘丰富区域
2. 对这些边缘丰富区域，提取相似区块分组，进行svd分解，选取大的特征值，重构得到去噪后的区块，并将得到的主成分作为词典典元进行保存
3. 利用以上的到的词典典元，计算相关系数（重构系数），对图像其余部分进行svd重构，实现去噪，若去噪效果不佳，则计算新的主成分，将其添加到词典中，更新词典，并以此重构
