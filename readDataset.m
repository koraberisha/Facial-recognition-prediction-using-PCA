function [vectorImages,susIds] = readDataset(folderName,rand,crop)

fullPath = pwd +"/"+ folderName;
vectorImages = zeros((112*92),400);
susIds = zeros(1,400);
count = 1;

for i=1:40
    k = 1;
    for n=1:10
        currentImage = fullPath + "/"+"s"+ i  + "/" + n + ".pgm";
        newImage = imread(currentImage);
        vectorImages(:,count) = newImage(:);
        susIds(count) = i;
        count = count+1;
        
    end

end



end


function vectorImages = generateRandomImages()

vectorImages = zeros((112*92),400);
sz = [10304, 400];
vectorImages = unifrnd(0,255,sz);


end
