% An example of how your function will be called, and what it should
% output.
% image_file_path is the absolute path to the image that you should
% process. This should be used to read in the file.
% image_file_name is just the name of the image. This should be written to
% the output file.
% output_file_path is the absolute path to the file where you should output
% the name of the file as well as the chocolates that you have detected.

% Uses Peter Corke's Machine Vision Toolbox for MATLAB.
function z3351846_MTRN4230_ASST1(image_file_path, image_file_name, ...
    output_file_path, program_folder)

im = imread(image_file_path);

chocolates = detect_chocolates(im);
chocolates
write_output_file(chocolates, image_file_name, output_file_path);

function c = detect_chocolates(im)
%% Background image
image = im;
I = rgb2gray(image);
% I = imadjust(I);
% index_pairs = [];

% Detect features
I_pts = detectSURFFeatures(I);
[I_features, I_validPts] = extractFeatures(I, I_pts);
figure;imshow(image);
% hold on; plot(I_pts.selectStrongest(50));

%% SUPER IMPORTANT VARIABLE
% empty_check = isempty(index_pairs);

%% Load reference image, and compute surf features
ref_img = imread('choc_up_ref.bmp');
ref_img_gray = rgb2gray(ref_img);
ref_pts = detectSURFFeatures(ref_img_gray);
[ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
index_pairs = matchFeatures(ref_features, I_features);
empty_check = isempty(index_pairs);

%% Ref image switching
% Switching to different ref image if different, meaning index_pairs is
% empty array.
% Failing both means choco not detected.
if empty_check
    ref_img = imread('choc_up_ref2.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    empty_check = isempty(index_pairs);
end;
if empty_check
    ref_img = imread('choc_down_ref.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    empty_check = isempty(index_pairs);
end;
if empty_check
    ref_img = imread('choc_down_ref2.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    empty_check = isempty(index_pairs);
end;
if empty_check
    ref_img = imread('choc_down_ref3.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    empty_check = isempty(index_pairs);
end;
if empty_check
    ref_img = imread('choc_down_ref4.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    empty_check = isempty(index_pairs);
end;
i = 1;

%% Showing points on ref image
% imshow(ref_img);
% hold on; plot(ref_pts.selectStrongest(50));

% %% Visual 25 SURF features
% % for showing the reference points, not really necessary
% figure;
% subplot(5,5,3); title('First 25 Features');
% for i=1:25
%     scale = ref_pts(i).Scale;
%     image = imcrop(ref_img,[ref_pts(i).Location-10*scale 20*scale 20*scale]);
%     subplot(5,5,i);
%     imshow(image);
%     hold on;
%     rectangle('Position',[5*scale 5*scale 10*scale 10*scale],'Curvature',1,'EdgeColor','g');
% end

%% SETTING BOUNDARIES 
% for determining reach and pick
r = 800;
circ_x = [0:1600];
circ_y= abs(220 + sqrt(r^2 - (circ_x - 800).^2));
    
% Safeguarding against blank space
if empty_check
    buffer_c(i,:) = zeros(1,8);
    buffer_c(i,7) = 1;
    corners = [];
    new_corners = [];
end;


%% While loop for multiple chocolates
while (empty_check ~= 1),
    %% Load reference image, and compute surf features
    ref_img = imread('choc_up_ref.bmp');
    ref_img_gray = rgb2gray(ref_img);
    ref_pts = detectSURFFeatures(ref_img_gray);
    [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
    index_pairs = matchFeatures(ref_features, I_features);
    face_indicator = 1;
    empty_check = isempty(index_pairs);
    % Switching to different ref image if different, meaning index_pairs is
    % empty array
    if (empty_check)
        ref_img = imread('choc_up_ref2.bmp');
        ref_img_gray = rgb2gray(ref_img);
        ref_pts = detectSURFFeatures(ref_img_gray);
        [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
        index_pairs = matchFeatures(ref_features, I_features);
        face_indicator = 1;
        empty_check = isempty(index_pairs);
    end;
    if isempty(index_pairs)
        ref_img = imread('choc_down_ref.bmp');
        ref_img_gray = rgb2gray(ref_img);
        ref_pts = detectSURFFeatures(ref_img_gray);
        [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
        index_pairs = matchFeatures(ref_features, I_features);
        face_indicator = 0;
        empty_check = isempty(index_pairs);
    end;
    if empty_check
        ref_img = imread('choc_down_ref2.bmp');
        ref_img_gray = rgb2gray(ref_img);
        ref_pts = detectSURFFeatures(ref_img_gray);
        [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
        index_pairs = matchFeatures(ref_features, I_features);
        face_indicator = 0;
        empty_check = isempty(index_pairs);
    end;
    if empty_check
        ref_img = imread('choc_down_ref3.bmp');
        ref_img_gray = rgb2gray(ref_img);
        ref_pts = detectSURFFeatures(ref_img_gray);
        [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
        index_pairs = matchFeatures(ref_features, I_features);
        face_indicator = 0;
        empty_check = isempty(index_pairs);
    end;
    if empty_check
        ref_img = imread('choc_down_ref4.bmp');
        ref_img_gray = rgb2gray(ref_img);
        ref_pts = detectSURFFeatures(ref_img_gray);
        [ref_features,  ref_validPts] = extractFeatures(ref_img_gray,  ref_pts);
        index_pairs = matchFeatures(ref_features, I_features);
        face_indicator = 0;
        empty_check = isempty(index_pairs);
    end;
    % Since index_pairs will be empty if there's no matches anymore
    % Will auto-break eventually at the end, however
    % ERROR on gte since it's expecting input.                                                                 
    % Error using vision.GeometricTransformEstimator/step
    % Empty on input 1 is not allowed.
    % Hence safeguard here.
    if empty_check
        break;
    end;

    ref_matched_pts = ref_validPts(index_pairs(:,1)).Location;
    I_matched_pts = I_validPts(index_pairs(:,2)).Location;
        figure, showMatchedFeatures(image, ref_img, I_matched_pts, ref_matched_pts, 'montage');
        title('Showing all matches');

    %% Define Geometric Transformation Objects
    %  Refining the image comparison, which works only if points  ref_matched_pts and I_matched_pts 3 <=
    % In case of false detection, will break.
    if((length(ref_matched_pts) < 3) && (length(I_matched_pts) < 3))
        if i == 1;
            buffer_c(i,:) = zeros(1,8);
            buffer_c(i,7) = 1;
            corners = [];
            new_corners = [];
        end;
        break;
    end;

    gte = vision.GeometricTransformEstimator;
    gte.Method = 'Random Sample Consensus (RANSAC)';

    [tform_matrix, inlierIdx] = step(gte, ref_matched_pts, I_matched_pts);

    ref_inlier_pts = ref_matched_pts(inlierIdx,:);
    I_inlier_pts = I_matched_pts(inlierIdx,:);

    % Draw the lines to matched points
    figure;showMatchedFeatures(image, ref_img, I_inlier_pts, ref_inlier_pts, 'montage');
    title('Showing match only with Inliers');

    %% Transform the corner points 
    % This will show where the object is located in the image

    tform = maketform('affine',double(tform_matrix));
    [width, height,~] = size(ref_img);
    corners = [0,0;height,0;height,width;0,width]
    new_corners = tformfwd(tform, corners(:,1),corners(:,2))
%     figure;imshow(image);
%     patch(new_corners(:,1),new_corners(:,2),[0 1 0],'FaceAlpha',0.5);
    
    %% Pickability indicator using equation of circle
    % Blocking out the strongest (detected) cluster
    drawRect = int32([new_corners(1,:) new_corners(2,:) new_corners(3,:) new_corners(4,:)]);
    shapeInserter = vision.ShapeInserter('Shape','Polygons','Fill', true, 'Opacity', 1);
    I = step(shapeInserter,I,drawRect);
    I_pts = detectSURFFeatures(I);
    [I_features, I_validPts] = extractFeatures(I, I_pts);
    figure;imshow(I);
    % patch(corners(:,1),corners(:,2),[0 1 0],'FaceAlpha',0.5);
    
    %% Storing the coordinates through x and y indicator updates
    x_indicator = 1600 - mean(new_corners(:,1));
    y_indicator = mean(new_corners(:,2));
    % Remember remember the flip of x axis
    local_delta_x = (new_corners(1,1) - new_corners(2,1));
    local_delta_y = (new_corners(1,2) - new_corners(2,2));
    ang_indicator = atan( local_delta_y / local_delta_x);
    if (local_delta_y > 0) && (local_delta_x > 0)
        ang_indicator = pi - ang_indicator;
    elseif (local_delta_y < 0) && (local_delta_x > 0)
        ang_indicator = -pi - ang_indicator;
    elseif ((local_delta_y < 0) && (local_delta_x < 0)) || ((local_delta_y > 0) && (local_delta_x < 0))
        ang_indicator = -ang_indicator;
    end;

    %% Pickable or not?
    % Index of circ_x == x_indicator
    foo = find(circ_x == round(x_indicator));
    if (y_indicator <= circ_y(foo))
        reach_indicator = 1;
        pick_indicator = 1;
    else
        reach_indicator = 0;
        pick_indicator = 0;
    end;
    %% Buffer update
    buffer_c(i,1) = x_indicator;
    buffer_c(i,2) = y_indicator;
    buffer_c(i,3) = ang_indicator;
    buffer_c(i,4) = 177;
    buffer_c(i,5) = 81;
    buffer_c(i,6) = face_indicator;
    buffer_c(i,7) = reach_indicator;
    buffer_c(i,8) = pick_indicator;

    i = i + 1;        
end;

%% For testing coords
hold on;
plot(circ_x,circ_y, 'b');

% th = linspace( 0, pi, 100);
% R = 800;  %or whatever radius you want
% circ_x = R*cos(th) + 800;
% circ_y = R*sin(th) + 220;
% plot(circ_x,circ_y, 'g'); axis equal;hold on;
% % preset_coord_x = [653.866584 800.500000 800.500000 ];
% % preset_coord_y = [564.589776 564.589776 564.589776];
% % plot(preset_coord_x,preset_coord_y,'*r');

if isempty(corners) ~= 1
% Coordinates of centroids
plot(corners(:,1), corners(:,2), 'r',new_corners(:,1), new_corners(:,2), 'b'); axis equal;
plot(buffer_c(:,1), buffer_c(:,2), '*g');
viscircles([1600 0], 10);
text(1400,0,'HERE BE ORIGINS')
text(corners(1,1), corners(1,2), '1');
text(corners(2,1), corners(2,2), '2');
text(corners(3,1), corners(3,2), '3');
text(corners(4,1), corners(4,2), '4');
text(new_corners(1,1), new_corners(1,2), '1');
text(new_corners(2,1), new_corners(2,2), '2');
text(new_corners(3,1), new_corners(3,2), '3');
text(new_corners(4,1), new_corners(4,2), '4');
end;



% You may store your results in matrix as shown below.
%    X    Y    Theta Width Height 1 = Top     1 = Reachable   1 = Pickable
%                                 0 = Bottom  0 = Unreachable 0 = Not pickable

c = buffer_c;

%  This is an example of how to write the results to file.
%  This will only work if you store your chocolates exactly as above.
%  Please ensure that you output your detected chocolates correctly. A
%  script will be made available so that you can run the comparison
%  yourselves, to test that it is working.
function write_output_file(chocolates, image_file_name, output_file_path)

fid = fopen(output_file_path, 'w');

fprintf(fid, 'image_file_name:\n');
fprintf(fid, '%s\n', image_file_name);
fprintf(fid, 'rectangles:\n');
fprintf(fid, ...
        [repmat('%f ', 1, size(chocolates, 2)), '\n'], chocolates');
    
% Please ensure that you close any files that you open. If you fail to do
% so, there may be a noticeable decrease in the speed of your processing.
fclose(fid);
% Three arguments:
%   full file path and name of input image
%   input image file name
%   input file image file name

% when called, this func will open image, process to calc info,
%   open results file and write results before closing borh files.


