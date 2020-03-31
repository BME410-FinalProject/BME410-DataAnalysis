addpath('D:\Documents\Word\College\Junior Year\Semester 2\BME 410\npy-matlab-master')
savepath

fn = 'D:\Documents\Word\College\Junior Year\Semester 2\BME 410\rawDataSample.bin';
fid = fopen(fn, 'r');
rawData = fread(fid, [385 Inf], '*int16');
fclose(fid);
chanMap = readNPY('channel_map.npy');
rawData = rawData(chanMap+1,:);
% figure; 
% imagesc(rawData(:,:));
figure;
hold on
for i = 1:5
    subplot(5, 1, i);
    plot(1705000:1800000, rawData(i, 1705000:1800000));
    set(gca, 'xticklabel', [])
end
hold off
