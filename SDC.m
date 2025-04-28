function [cl, icl] = SDC(rho, delta, aknn, nclust, aa)
n = length(rho);
r = rho.*delta;
k = size(aknn, 2);
cl = zeros(n, 1) - 1;
[~, ordr] = sort(r, "descend");
icl = [];

if nargin < 5
    aa = 0.4;
end

for i = 1:n
    if length(icl) == nclust
        break;
    end
    Ci = ordr(i);
    if cl(Ci) > 0
        continue;
    end
    icl = [icl; Ci];
    cl(Ci) = length(icl);
    queue = Ci;
    while ~isempty(queue)
        top = queue(1);  queue(1) = [];
        for j = 1:k
            xj = aknn(top,j);
            if cl(xj) > 0 || rho(xj) < aa*rho(Ci) || rho(xj) > (1+aa)*rho(Ci) %|| ~ismember(top, aknn(xj,:))
                continue;
            end
            queue = [queue; xj];
            cl(xj) = cl(Ci);
        end
    end
end
end
