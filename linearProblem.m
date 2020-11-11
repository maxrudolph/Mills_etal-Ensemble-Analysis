n = 3;
m = 3;
c = [5,4,3]';
b = [5,11,8]';
A = [2,3,1;4,1,2;3,4,2];

x = [zeros(n,1); b];
c = [c; zeros(m,1)];
A = [A, eye(m)];

setB = [n+1:n+m];
setN = [1:n];
B = A(:,setB);
N = A(:,setN);
xB = x(setB);
xN = x(setN);
cB = c(setB);
cN = c(setN);

xStarN = zeros(n,1);
xStarB = inv(B)*b;
zStarB = zeros(m,1);
zStarN = (inv(B)*N)'*cB - cN;
zetaStar = cB'*inv(B)*b;
zeta = zetaStar - zStarN'*xN;
xi = zetaStar + xStarB'*zB;

done = false;
if all(zStarN >= 0)
    done = true;
end
iteration = 0;
while ~done
    for iter = 1:n
        if zStarN(iter) < 0
            j = iter;
        end
    end
    j = setN(j);
    ej = zeros(n,1);
    ej(j) = 1;
    deltaXB = inv(B)*N*ej;
    
    options = zeros(m,1);
    for iter = 1:m
        options(iter) = deltaXB(iter)/xStarB(iter);
    end
    [t, i] = max(options);
    i = setB(i);
    t = 1/t;
    ei = zeros(m,1);
    ei(i) = 1;
    deltaZN = -(inv(B)*N)'*ei;
    s = zStarB(j)/deltaZN(j);
    
    xStarN(j) = t;
    xStarB = xStarB - t*deltaXB;
    zStarB(i) = s;
    zStarN = zStarN - s*deltaZN;
    
    for iter = 1:length(setB)
        if setB(iter) == i
            setB(iter) = [];
            break
        end
    end
    for iter = 1:length(setN)
        if setN(iter) == j
            setN(iter) = [];
            break
        end
    end
    setB = sort(setB, j);
    setN = sort(setN, i);
    
    if all(zStarN) >= 0
        done = true;
    end
    

iteration = iteration + 1

end

zeta
x(1:n)


 