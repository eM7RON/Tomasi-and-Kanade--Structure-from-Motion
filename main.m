% The layout of the code has been set up so you only need to run the top 5
% cells and should be able to generate any projection you like


% CLEAR VARIABLES AND TURN OFF IMAGE SIZE WARNING IF NEEDED

clear;
warning('off', 'Images:initSize:adjustingMag');
%%
% SELECT DATA FOR THE MODEL

pathway_string = "facevid%d\\facevid%d\\"; % Assumes unzipped file 
                                           % structure as provided
set_id = 1; % Select an image set (1 - 4)

% Only change one of the next two variables OR change neither of them
% to enable the use of all images to build the model.

input_image_ids = []; % File numbers of the images to be used to build the
                      % model. Leave empty to use all images or random
                      % sample.
                      
random_sample = 0; % If random sampling, set to size of random sample. 
                   % Else leave as 0 for no random sampling.

%%
% BUILD THE MODEL

[R, S, T, ids] = construct_model(...
                                 set_id,...
                                 pathway_string,... 
                                 input_image_ids,...
                                 random_sample...
                                 );

%%
% VIEW RESULTS

image_id = 1; % select image you would like to view. Be aware if 
              % using a subset this needs to be an index of the subset.
              
% Create 3D projection of image using model
              
r = R(image_id * 3 - 2: image_id * 3, :);
t = T(image_id * 3 - 2: image_id * 3, :);

X = r * S + t;

% get id of image to view

if ~isempty(ids)
    
    id = ids(image_id);
    
else
    
    id = image_id;
    
end

plot_sparse3Dface(X, set_id, id, ' (No processing)') % 3D projection of
                                                     % image
    
% comparison of original landmarks and projections in 2D
plot2D_landmarks_and_reprojections(...
                                   X,...
                                   set_id,...
                                   id,...
                                   pathway_string,...
                                   ' (No processing)'...
                                   );

%%
% TEST FOR BAS-RELIEF AMBIGUITY

if ambiguity(X)
        
    [R, S] = reconstruct_RS(R, S);

    r = R(image_id * 3 - 2: image_id * 3, :);
    t = T(image_id * 3 - 2: image_id * 3, :);
    
    X = r * S + t;
    
    % Retest ambiguity
    ambiguity(X);

    % 3D projection of image
    plot_sparse3Dface(X, set_id, id, ' (With processing)');

    % comparison of original landmarks and projections
    plot2D_landmarks_and_reprojections(...
                                       X,...
                                       set_id,...
                                       id,...
                                       pathway_string,...
                                       ' (With processing)'...
                                       );
end

%%
% Below are just my functions they do not need to be run.
% Some of the docstrings may refer to 'F' or 'P' which are the number
% of frames/images and number of points respectively

function [R, S, T, ids] = construct_model(...
                                          set_id,...
                                          pathway_string,...
                                          input_image_ids,...
                                          random_sample...
                                          )

    % This is a composite function which calls other functions.
    % You should be able to get an outline of the steps I have used.
    
    % Loads image data and constructs W
    [w, ids] = get_w(set_id, pathway_string, input_image_ids, random_sample);

    centroids = mean(w, 2); % Means of each row

    w = (w - centroids); % Centered for svd
    
    % Decomposition by SVD to yield estimate for R and S
    [Rh, Sh] = svd_w(w);

    % A system of linear equations is constructed from our R estimate
    [A, b] = construct_Ab(Rh);

    % Solves linear system to yield G
    G = get_G(A, b);

    % Eigen decomposition to yield Q. If this fails we will estimate Q
    % from SVD but this never seems to happen.
    Q = get_Q(G);
    
    % Yield true R and S using Q
    [R2, S] = solve_for_R2D_and_S(Rh, Sh, Q);
    
    % A new 3Fx3 R and 3Fx1 T are constructed
    [R, T] = construct_R_and_T(R2, centroids);
    
end

function [Rh, Sh] = svd_w(w)

    % Factorization of w by svd as described by Tomasi and Kanade 1992.

    [u, s, v] = svd(w);

    vt = v';

    s = sqrt(s(1: 3, 1: 3));

    Rh = u(:, 1: 3) * s;

    Sh = s * vt(1: 3, :);
    
end

function [R, S] = reconstruct_RS(R, S)

    % If Bas-relief ambiguity is detected, R will be corrected by
    % multiplying the 3rd column of each Ri by -1 and recalculating
    % the cross product to get 3rd row. The Z coordinates of S
    % will be multiplied by -1.
    
    for i = 1: size(R, 1) / 3
        
        R(i * 3 - 2: i * 3 - 1, 3) = R(i * 3 - 2: i * 3 - 1, 3) * - 1;
        
        R(i * 3, :) = cross(R(i * 3 - 2, :), R(i * 3 - 1, :)) * -1;
        
    end
    
    S(3, :) = S(3, :) * -1;
    
end

function out = ambiguity(X)

     % Detects bas-relief ambiguity by the assumption that the Z
     % coordinate of the nose tip is lower than the median Z 
     % coordinate.
     %
     % Args:
     %     X    - 2D array[float][float], projected 3D model of image, 
     %            shape[3, P]
     %
     % Returns:
     %      out - uint8, logical, 1 if ambiguity, 0 otherwise
     %          - disp, string to describe if ambiguity detected

     if X(3, 34) > median(X(3, :), 2)
        
        out = 1;
        
        disp("Bas ambiguity detected");
        
     else
        
        out = 0;
        
        disp("No Bas relief ambiguity detected");
        
    end
    
end

function [R, T] = construct_R_and_T(R2, centroids)

    % Constructs a 3F x 3 rotation matrix from a 2F x 3 rotation matrix
    % and a 3F x 1 translation vector from a 2F x 1 translation vector.
    % The new matrix is created where each pair of rows (x and y)
    % is preceeded by their cross product. The new translation matrix
    % is constructed from the centroids of the original input and
    % is padded with zeros between each x, y pair.
    
    n_rows = size(R2, 1) / 2 * 3;
    
    R = zeros(n_rows, 3);
    
    T = zeros(n_rows, 1);
    
    for i = 1: n_rows / 3
        
        R(i * 3 - 2: i * 3 - 1, :) = R2(i * 2 - 1: i * 2, :);
        
        R(i * 3, :) = cross(R2(i * 2, :), R2(i * 2 - 1, :));
        
        T(i * 3 - 2: i * 3 - 1, :) = centroids(i * 2 - 1: i * 2, :);
        
    end
    
end

function G = get_G(A, b)

    % Solve for X by linear least squares
    % G is then constructed from x.

    x = A \ b;
    
    G = [[x(1) x(2) x(3)];
         [x(2) x(4) x(5)];
         [x(3) x(5) x(6)]];
       
end
    

function [R2, S] = solve_for_R2D_and_S(Rh, Sh, Q)

    % Corrects Rh and Sh using Q.

    R2 = Rh * Q;
    
    S = Q \ Sh;
    
end

function Q = get_Q(G)
 
    Q = get_Q_eig(G);
    
    % If eigendecomposition not possible an estimate will be constructed
    % using svd.
    
    if ~Q
    
        disp("Deriving Q by eigendecomposition unsuccessful")
        
        Q = get_Q_svd(G);
        
        disp("Deriving estimate of Q by svd successful")
        
    else
        
    disp("Deriving Q by eigendecomposition successful")
        
    end

end

function Q = get_Q_svd(G)

    % This seems to yield a rough approximation for Q in that
    % projections look reasonable but lack some depth (Z magnitude)

    [u, s, v] = svd(G);
    
    Q = u(:, 1:3) * s(1: 3, 1: 3);

end

function Q = get_Q_eig(G)

    [u, l] = eig(G);
    
    if any(diag(l) < 0) % check for negative eigen values
        
        Q = 0;
        
    else   

        Q = u * sqrtm(l);
        
    end

end

function [A, b] = construct_Ab(Rh)

    % Construct A and b for solving to find G by linear least squares.
    % A will be a 3F x 6 array and b will be a 3F x 1 vector.

    F = size(Rh, 1) / 2;

    A = zeros(3 * F, 6);

    b = zeros(3 * F, 1);

    for i = 1: F
        
        x = Rh(2 * i - 1,:);

        y = Rh(2 * i,:);
        
        % Enforce constraints and orthogonality
        
        % i^T * i == 1 
        
        A(i * 3 - 2, :) = compute_row(x, x);
        
        b(i * 3 - 2, :) = 1;
        
        % j^T * j == 1 

        A(i * 3 - 1, :) = compute_row(y, y);

        b(i * 3 - 1, :) = 1;
        
        % i^T * j == 1 

        A(i * 3, :) = compute_row(x, y) * -1;

%       % b(i * 3, :) = 0 

    end

end

function row = compute_row(a, b)
    % Calculates a row for array A
    % As seen in:
    % A Sequential Factorization Method for Recovering Shape and Motion 
    % From Image Streams. Toshihiko Morita and Takeo Kanade, 
    % IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE,
    % VOL.  19,  NO.  8,  AUGUST  1997

row = [a(1) * b(1),...
       a(1) * b(2) + a(2) * b(1),...
       a(1) * b(3) + a(3) * b(1),...
       a(2) * b(2),... 
       a(2) * b(3) + a(3) * b(2),...
       a(3) * b(3)];
 
end







