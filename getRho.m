function [rho, aknn, aknn_dis] = getRho(data, numAnchor)

[n, m] = size(data);
generateAnchor = 1; % BKHK


clusts = AnchorGEN_noiter(data, numAnchor, generateAnchor);

kvalue = length(clusts{1});

clusts = cell2mat(clusts);
aknn = cell(n, 2);
for i = 1:size(clusts, 1)
    idx = clusts(i,:);
    sdist = EuDist2(data(idx,:), data(idx,:));
    for j = 1:kvalue
        aknn{idx(j), 1} = [aknn{idx(j),1} idx];
        aknn{idx(j), 2} = [aknn{idx(j),2} sdist(j,:)];
    end
end
aknn_dis = zeros(n, kvalue-1);
for i = 1:n
    [saknn, idx] = unique(aknn{i,1});
    item = aknn{i,2}(idx);
    [item, idx2] = sort(item);
    aknn{i,1} = saknn(idx2(2:kvalue));
    aknn_dis(i,:) = item(2:kvalue);
end
aknn(:,2) = [];
aknn = cell2mat(aknn);

% 高维数据中，欧式距离的距离度量会失效
% 来源于书《Introduction ot High-Dimensional Statistics-Chapman and
% Hall_CRC(2014)》, Figure 1.3
if m <= 100
    cv_thr = 0.4;
    m = 0;
    for i = 1:200
        m = 0.05 + m;
        rho = sum(aknn_dis(:, 1:end).^m, 2).^-1;
        cv = std(rho) / mean(rho);
        if cv > cv_thr
            break;
        end
    end
else
    rho = zeros(n, 1);
    for i = 1:n
        rho(aknn(i,:)) = rho(aknn(i,:)) + 1;
    end
end

end
