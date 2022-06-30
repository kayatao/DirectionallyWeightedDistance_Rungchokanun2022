clear all;

sourceFingerNetPath = 'D:\Research\FingerNet\FingerNet-master\output\20220120-095707\FVC2006_Db4FingerNet\FingerNetOriginal\';
destinationTextPath = 'D:\Research\FingerNet\FingerNet-master\output\20220120-095707\FVC2006_Db4FingerNet\OF\';
files = dir(fullfile(sourceFingerNetPath, '*.mat'));

for i = 1 : size(files,1)
    load(fullfile(sourceFingerNetPath, files(i).name))
    fileIDF = fopen(destinationTextPath + string(strrep(files(i).name,'.mat','.txt')),'w');
    for j = 1 : size(orientation,2)
        fprintf(fileIDF,'%f ',orientation(:,j,:));
        fprintf(fileIDF,'\n');
    end
    fclose(fileIDF);
    clear orientation orientation_distribution_map
end