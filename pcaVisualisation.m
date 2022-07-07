
function pcaVisualisation()
A = readDataset("att_faces",1);


[V,L,mu] = cw_pca(A);
D = dct2(A),10304,400;

imagesc(reshape(mu,112,92));

pcaBasisDisplayer(5,V);
pcaDisplayEigenenergy(V,L);

imageTest = A(1:(112*92),1);

basisTest = V(1:(112*92),1);

projectScatter2(A,V,mu);

newV = trancateBasisProjection(A,V,20,mu);

end

function pcaBasisDisplayer(noOf,V)

imagesc([reshape(V(:,1),112,92),reshape(V(1:(112*92),2),112,92),reshape(V(1:(112*92),3),112,92),reshape(V(1:(112*92),4),112,92),reshape(V(1:(112*92),5),112,92)]);

end

function pcaDisplayEigenenergy(V,L)

energy = 0;
arr = zeros(1,400);
for n=1:400
    energy = energy + L(n);
    for i=1:5
        current = V(:,i);
        arr(1,n) = energy;
        
    end
end

plot(arr);

end


function [p,xn] = projectVecs(x,V,mu,k)
    
Vd = V.';
% Vd = round(Vd*1000);
% Vd = Vd / 1000;
p = (Vd*(x-mu));
xn = (V*p+mu);

    
end


function vals = systematicallyAlter(x,V,mu,L,index)
ai = 11
aj = 0
v = 0
newArr = V(1:10304,index);
variation = newArr * mu+(v*sqrt(L(index))) ;

vals = unprojectVecs(x,variation,mu);

end


function projectScatter2(A,V,mu)

p1 = [];
p2 = [];
for n=1:400
[p,val] = projectVecs(A(:,n),V,mu);
p1(n) = p(2);
p2(n) = p(5);
    
    
end
c = linspace(1,10,length(p1));
scatter(p1,p2,50,c,'filled');
colormap parula;



end

function imrec = quantisationMethod(A)

im = A;
[h, w] = size(im);
len = numel(im);

x1 = im(:, 1:2:end);
x2 = im(:, 2:2:end);
X = [x1(:) x2(:)]';

Q = 64; % Quantisation factor (for each component of the pair, so overall factor = Q*Q)

X = X .* 255;
Xq = round(X ./ Q);
% To transmit:
% Quantization constant, Q
% Codes, Xq (nearest centres)

% Reconstruct the image from codes
imrec = zeros(size(im));
imrec(:, 1:2:end) = reshape(Xq(1, :) * Q, h, w/2);
imrec(:, 2:2:end) = reshape(Xq(2, :) * Q, h, w/2);
imrec = imrec ./ 255;


end

function ratio = calculateCompression(origIm,V,p,mu)

origIm = whos("origIm");
V = whos("V");
p = whos("p");
mu = whos("mu");


A = origIm.bytes;
B = (V.bytes*p.bytes)+mu.bytes;

ratio = B/A;


end

function createCompressionPlot(x,mu,V)
compressions = zeros(1,400);

for n=1:40
    
    m = n*10;
     curr = V(:,1:m);
    [p,xn] = projectVecs(x,curr,mu,1);
    
    ratio1 = calculateCompression(x,V,p,mu)
    compressions(m) = ratio1;
end

vals = linspace(1,400,400);
plot(vals,compressions);

end


function newV = trancateBasisProjection(x,V,k,mu)


newV = 0;
[p,xn] = projectVecs(x(:,1),V(:,1:10),mu,1);
a1 = reshape(xn,112,92);
x1 = xn;
[p,xn] = projectVecs(x(:,1),V(:,1:50),mu,1);
a2 = reshape(xn,112,92);
x2 = xn;
[p,xn] = projectVecs(x(:,1),V(:,1:250),mu,1);
a3 = reshape(xn,112,92);
x3 = xn;
[p,xn] = projectVecs(x(:,1),V(:,1:350),mu,1);
a4 = (reshape(xn,112,92));
x4 = xn;
b = [a1,a2,a3,a4];
xnArr = [x1,x2,x3,x4];


xa = zeros(4);

for n=1:4
   xa(n) = sqrt(immse(x(:,1), xnArr(:,n)));
    
end



plot(xa)

% endval = 10304/4;
% [p,xn1] = projectVecs(x(1:endval,1),V(1:endval,1:400),mu(1:endval),1);
% 
% [p,xn2] = projectVecs(x(endval:endval*2,1),V(endval:endval*2,1:400),mu(endval:endval*2),1);
% 
% [p,xn3] = projectVecs(x(endval*2:endval*3,1),V(endval*2:endval*3,1:400),mu(endval*2:endval*3),1);
% 
% [p,xn4] = projectVecs(x(endval*3:end,1),V(endval*3:end,1:400),mu(endval*3:end),1);
% 
% 
% 
% newVal = zeros(1,10304);
% newVal(1:endval) = xn1;
% newVal(endval:endval*2) = xn2;
% newVal(endval*2:endval*3) = xn3;
% newVal(endval*3:end) = xn3;
% size(newVal)
% 
% imageDiscect = reshape(newVal,112,92);
% imdisp(imageDiscect);


end


