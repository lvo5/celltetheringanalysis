%% Read the video file
a=VideoReader('EC6.avi'); %Name of video here.
%Write file into a folder
mkdir('EC6'); %Write the folder name as your video name 
filepath = '/Users/lvo5/Desktop/code for tethering/EC6/';
for img = 1:a.NumberOfFrames;
    framename= strcat('frame',num2str(img),'.tiff');
    filename = strcat(filepath, framename);
    b = read(a, img);
    b = rgb2gray(b);
    imwrite(b,filename);
end
%% Gage the binary condition to determine threshold parameter
% Gage cropping:  This step is important because it is useful if you binary image
% shown above have alot of background noises. This will help the analysis
% goes faster!
 
picture = strcat(filepath,'frame1.tiff');
grayImage = imread(picture);
subplot(2, 2, 2); 
histogram(grayImage, 0:256);
grid on;
title('Histogram of original image', 'FontSize', 13, 'Interpreter', 'None');
xlabel('Gray Level', 'FontSize', 13);
ylabel('Pixel Count', 'FontSize', 13);
xlim([0 255]); 

% Scale x axis manually.
% Threshold and binarize the image

% You can adjust until you are happy with the cell binary
% Display the image.

binaryImage = grayImage > 95; % Adjust your binary condition.

%Here you will have to pick the pixel coordinate of your cells. You would
%have to pick the cell that are rotating from a black and white picture.
% This can be difficult; however, it will be worth it at the end. This is
% the only hard part of the code. 

%The coordinate system should be in form of 
%Remove the % for the line below.

%binaryImage =  imcrop(binaryImage,[[1680, 609, 1710-1680, 640-610]]);

% the coordinate should be in form of [x1 y1 delx1 dely1]

subplot(2, 2, 3);
imshow(binaryImage, []);
axis on;
title('Binary Image', 'FontSize', 13, 'Interpreter', 'None');
labeledImage = bwlabel(binaryImage);
%% convert to binary images and crop the desired cell.
% Here you would have to put all of your coordinates in form of
%[[x1 y1 delx1 dely1]; [x2 y2 delx2 dely2], ...]
% For this example, you have two cells rotating so i'm gonna put two.
% [1680, 609, 1710-1680, 640-610]

%Accounting for uneven binary threshold: here you can put the cell number
%that does not cannot be caught with the threshold of the rest of the cell

% For this example, cell #3 at 101 threshold cannot be seen. But can be
% seen clearly at 95 threshold.

coordinates_box = [[1180, 510, 1240-1180, 570-510]; [570, 290, 630-570, 330-290]; [1680, 609, 1710-1680, 640-610]];
matrix_size = size(coordinates_box); 

%Here put in the number of cell that is different:
cell_diff = [3]; %in form of [1, 2, 3....]
threshold_diff = [95]; %threshold matrix for exceptions in form of [101, 95, ...]

%This example means threshold 95 for cell 3. 


for i2 = 1:matrix_size(1);
    format = 'EC6_binary_%d'; 
    folder_binary = sprintf(format, i2);
    mkdir(folder_binary);
    filepath = '/Users/lvo5/Desktop/code for tethering/EC6/';
    
    if ismember(i2, cell_diff);
        index_dif = find(cell_diff == i2);
        
        for i = 1:a.NumberOfFrames;
    %read file format
            format = 'frame%d.tiff';
            framename = sprintf(format, i);
            filename = strcat(filepath, framename);
            grayImage = imread(filename);
            binaryImage = grayImage > threshold_diff(index_dif); %change threshold based on previous section
            binaryImage =  imcrop(binaryImage, coordinates_box(i2,:));
    
            formatsave = 'frame%d_binary.tiff';
            filename_save = sprintf(formatsave, i);
        %Put into seperate folders of binary images
            foldername = '/Users/lvo5/Desktop/code for tethering/';
            foldername = strcat(foldername,  folder_binary);
            fullFileName = fullfile(foldername, filename_save);
            imwrite(binaryImage,fullFileName);
        end 
    else
        for i = 1:a.NumberOfFrames;
    %read file format
            format = 'frame%d.tiff';
            framename = sprintf(format, i);
            filename = strcat(filepath, framename);
            grayImage = imread(filename);
            binaryImage = grayImage > 101; %change threshold based on previous section
            binaryImage =  imcrop(binaryImage, coordinates_box(i2,:));
    
            formatsave = 'frame%d_binary.tiff';
            filename_save = sprintf(formatsave, i);
        %Put into seperate folders of binary images
            foldername = '/Users/lvo5/Desktop/code for tethering/';
            foldername = strcat(foldername,  folder_binary);
            fullFileName = fullfile(foldername, filename_save);
            imwrite(binaryImage,fullFileName);
        end 
    end
end
%% Checkpoint: make sure that rotating cell are captures, we can
% make a movie

for i2 = 1:matrix_size(1);
   %Open video file with names
    videoname = strcat('WT_sample_', num2str(i2)); % Change video name to anything
    videoname = strcat(videoname, '.MP4'); %change filetype if you want
    movie_obj = VideoWriter(videoname); %Write any name here
    open(movie_obj);
    
    %Make video
    format = 'EC6_binary_%d'; %change name here
    folder_binary = sprintf(format, i2);
    filepath = '/Users/lvo5/Desktop/code for tethering/';
    hi = strcat(filepath,  folder_binary);
    filepath = strcat(hi,  '/');
    for K = 1 : a.NumberOfFrames;
        frame_template = 'frame%d_binary.tiff';
        filename = sprintf(frame_template, K);
        filename = strcat(filepath, filename);
        this_image = imread(filename);
        this_image = im2double(this_image);
        writeVideo(movie_obj, this_image);
    end
end

close(movie_obj);

%% Making matrix for all objects detected making data: Here we collect 
% the raw coordinate data. Sometime, each frame will be contaminated with
% raw pixels from back. To eliminate this, we can set the pixel threshold.
% The default is 12, but you can change it. Ideally, you want one object in
% all frames. Pixel threshold is just the major axis of the ellispes, which
% is what these cells are modeled by.
% Hierachy (top to bottom)  -> empty, cell#, frame#, variables)

data = {};
pixel_threshold = 12; % Change pixel threshold for major axis
for i = 1:matrix_size(1);
    dummy_matrix = [];
    cell_name = 'EC6_binary_';
    cell_name = strcat(cell_name, num2str(i));
    for i2 = 1:a.NumberOfFrames;
        frame_template = 'frame%d_binary.tiff';
        filename = sprintf(frame_template, i2);
        
        filepath = '/Users/lvo5/Desktop/code for tethering/';
        filepath = strcat(filepath, cell_name);
        filepath = strcat(filepath, '/');
        
        A=imread(strcat(filepath, filename));
        A = imcomplement(A);
        label = bwconncomp(A, 4);
        measurements = regionprops(label, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
        for i3 = 1: length(measurements);
            if measurements(i3).MajorAxisLength > pixel_threshold;
             measurements_isolated = measurements(i3);
            end
        end
    % Adding variable 
    dummy_matrix(i2,1) = measurements_isolated.Centroid(1); 
    dummy_matrix(i2,2) = measurements_isolated.Centroid(2);
    dummy_matrix(i2,3) = measurements_isolated.MajorAxisLength; 
    dummy_matrix(i2,4) = measurements_isolated.MinorAxisLength;
    dummy_matrix(i2,9) = i2* a.Duration/a.NumberOfFrames; 
    end
    data{i} = dummy_matrix; 
end

% Setting another variable for raw:
data_raw = data;
%% Fixing outliers from dataset for better fixation from least square regression
for i = 1:matrix_size(1);
    %Fixing the outliers through statistics
    B_x = filloutliers(data_raw{i}(:,1),'center','quartiles');
    B_y = filloutliers(data_raw{i}(:,2),'center','quartiles');
    
    data_raw{i}(:,1) =  B_x; 
    data_raw{i}(:,2) =  B_y;
end

%% Using the least square fit method,find the approximated center of the
% centroid coordinates. I used the code that is already made for me
% https://www.mathworks.com/matlabcentral/fileexchange/22643-circle-fit-pratt-method

%--------------------------------------------------------------------------
%  
%     Circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this fit does not use built-in matrix functions (except "mean"),
%           so it can be easily programmed in any programming language
%
%--------------------------------------------------------------------------
 
for i = 1:length(data_raw); %Here, we fit the center of the dataset
   circlefit = CircleFitByPratt(data_raw{i}(:,[1,2]));
   for i2 = 1:length(data_raw{i});
       data{i}(i2, 6) = circlefit(1) ; %fitted x center
       data{i}(i2,7) = circlefit(2); %fitted y center
       data{i}(i2, 8) = circlefit(3); %radius of rotation
   end
end

%Apply the moving average smoothing procedure.
for i = 1:matrix_size(1);
    data{i}(:,1) = movmean(data{i}(:,1),2);
    data{i}(:,2) = movmean(data{i}(:,2),2);
end

%% Calculate angle of rotation: we can use the following formula:
% tan_inv((centroidy-centery)/(centroidx-centerx));

for i = 1:length(data);
    diffx = data{i}(:,1) - data{i}(:,6);
    diffy = data{i}(:,2) - data{i}(:,7);
    angle = atan2(diffy, diffx); 
    %Column #5 is the true angle of rotation
    winddir = angle;
    data{i}(:,5) = winddir * 180/pi; 
    data{i}(:,10) = unwrap(winddir)* 180/pi; %Unwrap  
end 
%% Checkpoint: Here you want to remove any bad tracks!
%Here, we are going to graph all of the track selected and the fit. This
%Step is to determine whether the dataset for each cell is good or not. We 
% can discard the bad one on the next section.

% For this example, the second cell is bad. So, you have to remove in the
% next section
close all
degree = [0:1:360];
for i = 1:matrix_size(1);
    figure()
    plot(data{i}(:,1),data{i}(:,2));
    hold on
    plot(data{i}(1,8)*cosd(degree)+data{i}(1,6),data{i}(1,8)*sind(degree)+data{i}(1,7),'red');
    scatter(data{i}(1,6), data{i}(1,7), 'k');
    title(strcat('EC6cell#', num2str(i)));
    
    axis equal
    hold off
end 
%% Discarding the bad cells - OPTIONAL: if you see that one of the scatter
% dataset is bad, you can remove them (which I suggest).

idcell = [1,2,3]; %Here you put the CELL NUMBER ID YOU WANT TO KEEP!
% put in increasing order.

data = data(idcell);

%% Basic statistics: calculate revolution per minutes, angular velocity, stop frequencies, reversal frequency per minutes. 
for i = 1:length(data);
    % Calculate velocity profile over time:
    displacement = diff(data{i}(:,10));
    dt = diff(data{i}(:,9));
    angular_velocity = displacement./dt; 
    angular_velocity1 = vertcat([NaN], angular_velocity);
    data{i}(:,11) = angular_velocity1;
    
    %Angular acceleration  
    dv= diff(angular_velocity);
    dt1 = dt(1,:);
    angular_acceleration = dv./dt1;
    angular_acceleration = vertcat([NaN; NaN], angular_acceleration);
    data{i}(:,13) = angular_acceleration;
    
    %Revolution per minutes (speed)
    speed = abs(angular_velocity) * pi/180;
    speed = speed/(2*pi);
    speed = vertcat([NaN], speed);
    data{i}(:,12) = speed;
    
    pauses = [];
    pause_threshold = 0.2; %speed limit 0.2 rad/sec is a good threshold
    %Create a for matrix for pause detection. 
    for i2 = 1:length(speed);
        if speed(i2) <= pause_threshold;
            pauses(i2) = 1;
        else
            pauses(i2) = 0;
        end
    end
    data{i}(:,14) = pauses;
    
    %Direction of rotation 
    direction = [];
    for i2 = 1:length(angular_velocity1);
        if angular_velocity1(i2) < 0;
            direction(i2)= -1;
        else
            direction(i2) = 1;
        end
    end
    data{i}(:,15) = direction;
end

% Creating a summary matrix: for all cells.
summary = []; %cell#, averagespeed(Hz),averagereversalfrequency,  averagepausefrequency,...
                %averagepauseduration
for i = 1:length(data);
    %Cell numbers
    summary(i, 1) = i; 
    
    %Calculate average speeds (Hz)
    avg_speed = nanmean(data{i}(:,12));
    summary(i,2) = avg_speed;
    
    %Calculate reversal frequency: Count sign changes and divide by total time
    frequency = data{i}(:,15);
    change_count = 0;
    for i2 = 1:length(frequency);
        if i2 == length(frequency);
            break
        end
        
        if frequency(i2) * frequency(i2 + 1) < 0;
            change_count = change_count + 1;
        end
    end
    summary(i,3) = change_count./a.Duration;
    
    %Calculate average pausing frequency and average pause time
    pause = data{i}(:,14);
    [r,s] = runlength(pause,numel(pause));
    pause_matrix = r(find(s==1));
    
    %Pause frequency:
    summary(i,4) = length(pause_matrix)./a.Duration;
    
    %Pause duration:
    summary(i,5) = mean(pause_matrix.*dt(1));
end 


%% Save raw dataset in excel file

filename = 'EC6_rawdata.xlsx' ; %Name anything you want. 
filepath = '/Users/lvo5/Desktop/code for tethering/';
for i= 1:length(idcell);
    size1 = size(data{i});
    header = {'CentroidX', 'CentroidY', 'MajorAxisLength', 'MinorAxisLength','Angle of rotation', 'Fitted_center_x', 'Fitted_center_y', 'Radius (pixel)' , 'Time (sec)', 'Unwrapped angle','Angular velocity (deg/sec)',...
        'Speed (Hz)', 'Angular acceleration (deg/sec^2)', 'Pause events', 'Direction'};
    c = cell(size1(1)+1, size1(2)); 
    hi = size(c);
    
    c(1,:) = header;
    c(2:hi(1),:) = num2cell(data{i});
    
    substring = strcat('cell#', num2str(idcell(i)));
    writecell(c, strcat(filepath, filename) ,'sheet' , substring)
end

%Summary dataset
size1 = size(summary);
header = {'ID', 'Average Hz', 'Avg_reversalfreq', 'Avg_pausefreq', 'Avg_pausetime' };
c = cell(size1(1)+1, size1(2));
hi = size(c);
c(1,:) = header;
c(2:hi(1),:) = num2cell(summary);
writecell(c, strcat(filepath, filename) ,'sheet' , 'Summary')

%% Plotting example 
for i = 1:length(data);
figure()
subplot(2, 1, 1)
plot(data{i}(:,9) , data{i}(:,10), 'blue')
title(strcat('Cell#', num2str(i)))
ylabel('Degree')
hold on 
%Scattering pauses event 
index = find(data{i}(:,14) == 1);
time = data{i}(:,9);
pauset = time(index);
hi = data{i}(:,10);
scatter(pauset , hi(index), 'red', '*')
hold off

subplot(2,1,2)
plot(data{i}(:,9) , data{i}(:,11), 'black')
hold on
hi = data{i}(:,11);
scatter(pauset , hi(index), 'red', '*')
ylabel('Angular velocity (degree/sec)')
xlabel('Time (sec)')
%Scattering pauses event 
hold off
end
