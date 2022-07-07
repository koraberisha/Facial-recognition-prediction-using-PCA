function eigenfaceRecognition()
ids = zeros(1,400);

count = 1;
nval = randi([1,400],1);


noOf = 0;
[AA,sus] = readDataset("att_faces",1);

for nval=1:400
A=AA;
value = A(:,nval);
A(:,nval) = [];
[V,L,mu] = cw_pca(A);
[projCoefficients,featureVec] = projectVecs(A,V,mu);
[suspectCoefficients,suspectFeature] = projectVecs(value,V,mu);
% projectScatter2(projCoefficients,suspectCoefficients);

distFeature = featureVec';
distSuspect = suspectFeature';

[D,I] = pdist2(distFeature,distSuspect,"euclidean","Smallest",3);
if (sus(nval) == sus(I(1)))
    noOf = noOf + 1;

end

end

noOf 
imdisp([reshape(value,112,92),reshape(A(:,(I(1))),112,92),reshape(A(:,(I(2))),112,92)])

end

function vVal = findDist(u,v)
newLow = 10000;
nVal = 0;
vVal = [];

for n=1:399
    CosTheta = max(min(dot(u(:,n),v)/(norm(u(:,n))*norm(v)),1),-1);
    ThetaInDegrees = real(acosd(CosTheta));
    vVal(n) = ThetaInDegrees;
end

end

function projectScatter2(fet,sus)

p1 = [];
p2 = [];
for n=1:399
p1(n) = fet(2,n);
p2(n) = fet(3,n);
        
end

p1(400) = sus(2);
p2(400) = sus(3);

c = linspace(1,100,length(p1));
scatter(p1,p2,30,c,'filled');
colormap jet;


end

function u = eigenfaceGeneration(V,L)

u = zeros(10304,399);
for n=1:399   
    u(:,n) = V(:,n)*L(n);
end

end

function [p,xn] = projectVecs(x,V,mu)
    
Vd = V.';
p = (Vd*(x-mu));
xn = V*p+mu;
    
end