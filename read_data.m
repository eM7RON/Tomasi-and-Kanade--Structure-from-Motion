


function w = read_data(n_file, n_images)

% A function to read in the image landmark data
% Please be aware you may have to adjust the path based on your own file
% structure - this function is designed to use the unzipped file structure
%
%    Args:
%        n_file   - uint8, which set of images is to be loaded 1 - 4
%        n_images - uint8, number of images in set

proto_name = sprintf("facevid%d\\facevid%d\\", n_file, n_file);

    file_extension = ".pts";

    x = [];

    y = [];
    
    i = 1;

    while 1

        if i < 10

            file_number = strcat("000", num2str(i));

        elseif i < 100

            file_number = strcat("00", num2str(i));

        else

            file_number = strcat("0", num2str(i));

        end

        file_name = strcat(proto_name, file_number, file_extension);
        
        disp(file_name);
        
        try

            intermediate = read_landmarks(file_name);
            
            x = [x; intermediate(:, 1)];
        
            y = [y; intermediate(:, 2)];
            
        catch
            
            break
            
        end
        
        i = i + 1;
        
    w = [x; y];
             
    end
    
   
   
    
        

    
