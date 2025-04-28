function cl = FDPAN(data, nclust, numAnchor, k)

data = normalize(data, "range");
[n, m] = size(data);

[rho, aknn, aknn_dis] = getRho(data, numAnchor);

if n <= 20000 && m <= 2
    [aknn, aknn_dis] = knnsearch(data, data, "K", k+1);
    aknn(:, 1) = [];  aknn_dis(:, 1) = [];
else
    aknn = aknn(:, 1:k);  aknn_dis = aknn_dis(:, 1:k);
end

nneigh = zeros(n, 1);
delta = zeros(n, 1);
for i = 1:n
    for j = 1:k
        jnn = aknn(i,j);
        if rho(jnn) > rho(i) && ismember(i, aknn(jnn,:))
            nneigh(i) = jnn;
            delta(i) = aknn_dis(i,j);
            break;
        end
    end
end

[~, ordrho] = sort(rho, "descend");
icl = [];
for i = 1:n
    if nneigh(ordrho(i)) == 0
        icl = [icl ordrho(i)];
    end
end

numC = length(icl);
if numC < nclust
    cl = 0;
    return;
end

level = numAnchor;
while numC > 20000
    [~, aknn_icl, aknn_dicl] = getRho(data(icl,:), level);
    while size(aknn_icl, 2) < 5
        level = level - 1;
        [~, aknn_icl, aknn_dicl] = getRho(data(icl,:), level);
    end
    for i = 1:numC
        for j = 1:size(aknn_icl, 2)
            jnn = aknn_icl(i,j);
            if rho(icl(jnn)) > rho(icl(i)) && ismember(i, aknn_icl(jnn,:))
                nneigh(icl(i)) = icl(jnn);
                delta(icl(i)) = aknn_dicl(i,j);
                break;
            end
        end
    end
    icl = icl(nneigh(icl)==0);
    numC = length(icl);
end
clear aknn_dis 

distC = EuDist2(data(icl,:), data(icl,:));
[~, ordC] = sort(rho(icl), "descend");
for i = 2:numC
    delta(icl(ordC(i))) = inf;
    for j = 1:i-1
        if distC(ordC(i), ordC(j)) < delta(icl(ordC(i)))
            delta(icl(ordC(i))) = distC(ordC(i), ordC(j));
            nneigh(icl(ordC(i))) = icl(ordC(j));
        end
    end
end
delta(icl(ordC(1))) = max(delta);

% figure; plot(data(:,1), data(:,2), 'b.');
[cl, ~] = SDC(rho, delta, aknn, nclust);
for i = 1:n
    if cl(ordrho(i)) < 0
        cl(ordrho(i)) = cl(nneigh(ordrho(i)));
    end
end
end