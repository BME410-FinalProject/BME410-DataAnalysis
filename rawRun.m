fn = 'D:\Documents\Word\College\Junior Year\Semester 2\BME 410\rawDataSample.bin';
fid = fopen(fn, 'r');
dat = fread(fid, [385 Inf], '*int16');
fclose(fid);
chanMap = readNPY('channel_map.npy');
dat = dat(chanMap+1,:);
figure; imagesc(dat(:,1:30000))