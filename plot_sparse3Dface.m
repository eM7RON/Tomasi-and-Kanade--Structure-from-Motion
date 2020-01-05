function plot_sparse3Dface(Shape, set_id, image_id, title_arg)
% Plot the sparse 3D face reconstruction in an intuitive way, exploiting
% the specific anatomical structure of the face landmarks: Based on the
% iBug68 landmark indexing, plot different sub-sets of landmarks as lines
% (e.g.~left eyebrow as an open polygonal line, left eye as a closed
% polygonal line, etc.). 
% This function uses Red colour for the Right eye, green for the left eye
% and blue for everything else. 
%
% (INPUT):
% Shape: 3 x Npoints array with the 3D coordinates of the sparse 3D facial
% shape to visualise. The points should follow the iBug68 markup system.
% Each column includes the XYZ coordinates of the corresponding landmark
% point
%
% Note: please feel free to modify this function according to your needs

% set of landmarks to be visualised as lines:
lines_lands_vis = {1:17, 18:22, 23:27, 28:31, 32:36, [37:42 37], [43:48 43], [49:60 49], [61:68 61]};

% specify the properties of the line to draw per segment, so that the right
% eye (indices 37:42 in iBug markup) is RED and the left eye (indices 43:48
% in iBug markup) is GREEN
line_properties = {'.-b', '.-b', '.-b', '.-b', '.-b', '.-r', '.-g', '.-b', '.-b'};

figure;
for i_line = 1:length(lines_lands_vis)    
    plot3(Shape(1,(lines_lands_vis{i_line})), Shape(2,(lines_lands_vis{i_line})), Shape(3,(lines_lands_vis{i_line})), line_properties{i_line}, 'LineWidth', 1.5);
    hold on;
end
hold off;
axis equal;

grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');

title_string = sprintf('3D projection for image %d of set %d', image_id, set_id);

title_string = strcat(title_string, title_arg);

title(title_string);

cameratoolbar('Show'); cameratoolbar('SetMode','orbit'); cameratoolbar('SetCoordSys','none');
end