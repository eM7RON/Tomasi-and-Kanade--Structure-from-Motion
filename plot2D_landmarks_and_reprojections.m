
function plot2D_landmarks_and_reprojections(...
                                            X,...
                                            set_id,... 
                                            image_id,...
                                            pathway_string,...
                                            title_arg...
                                            )

% outputs a figure of landmarks plotted on corresponding image
%
% Args:
%    set_id   - uint8, the set the image belongs to
%    image_id - uint8, the number of the image in the set
%    R        - array[float][float], rotation matrix, shape[3 * F, 3]
%    S        - array[float][float], shape matrix, shape[3, P]
%    T        - array[float], translation vector, shape[3 * F, 1]
%    (where F and P are the total number of images and total number of 
%     landmarks used to create R, S and T respectively).
%
% Output:
%    figure   - outputs a figure of landmarks plotted on corresponding
%               image

proto_name = sprintf(pathway_string, set_id, set_id);

pic_extension = ".jpg";

lm_extension = ".pts";

if image_id > 99

    lm_id = strcat("0", num2str(image_id));
    im_id = strcat("000", num2str(image_id));

elseif image_id > 9
    
    lm_id = strcat("00", num2str(image_id));
    im_id = strcat("0000", num2str(image_id));

else

    lm_id = strcat("000", num2str(image_id));
    im_id = strcat("00000", num2str(image_id));

end

% get image
image_name = strcat(proto_name, im_id, pic_extension);

image = imread(image_name);

% get landmarks
lm_name = strcat(proto_name, lm_id, lm_extension);

lm = read_landmarks(lm_name);

figure;
 
showMatchedFeatures(image, image, lm, X(1: 2, :)', 'blend', 'PlotOptions', {'ro','g+','b-'});

title_string = sprintf('Landmarks (o) and Projections (+) for image %d of set %d', image_id, set_id);

title_string = strcat(title_string, title_arg);

title(title_string);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.5, 0.5, 0.5]);

end

