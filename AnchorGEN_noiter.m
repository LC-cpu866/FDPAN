% AnchorGEN()的非递归版本
function clusts = AnchorGEN_noiter(X, numAnchor, selAnchor)

[n, ~] = size(X);
overlap = ceil(n / 2^(numAnchor)); % 重叠度设置为平衡划分中，点数或点数的一半

clusts = cell(2^numAnchor, 1);
clusts{1} = (1:n);
for level = 1:numAnchor
    numc = 2^(level-1);
    for i = 1:numc
        idx0 = clusts{i};
        y = kmeans_two(X(idx0,:), overlap);
        clusts{i} = idx0(y(:,1) > 0);
        clusts{numc + i} = idx0(y(:,2) > 0);
    end
end

% clusts = clusts';
end

% 二分类kmeans
function F = kmeans_two(X, overlap)

class_num = 2;
n = size(X, 1);
ratio = 0.5;

StartInd = randsrc(n, 1, 1:class_num);
InitF = TransformL(StartInd, class_num);

a = floor((n+1) * ratio) + overlap;

F = InitF;
for iter = 1:10
    C = X'*F./sum(F + eps);
    F = zeros(n, class_num);
    Q = EuDist2(X, C', 0);
    % Q = py.scipy.spatial.distance.cdist(py.numpy.array(X), py.numpy.array(C'), 'euclidean');
    q = Q(:,1) - Q(:,2);
    [~, idx] = sort(q);

    cp = a;
    if cp < 1
        cp = 1;
    elseif cp > n-1
        cp = n-1;
    end

    F(idx(1:cp), 1) = 1;
    F(idx(n-cp+1:n), 2) = 1;
    % 算法收敛时，跳出循环
    if norm(InitF - F, "fro") < 1e-6
        break;
    end
    InitF = F;
end
end