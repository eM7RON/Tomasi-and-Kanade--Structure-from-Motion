
function [w, ids] = get_w(...
                          set_id,...
                          pathway_string,... 
                          ids,...
                          random_sample...
                          )

    % A function to read in the image landmark data.
    % Please be aware it relies on the file structure created by unzipping
    % the files provided.
    %
    % Args:
    %     n_file          - uint8, which set of images is to be loaded 1 - 4
    %     pathway_string  - string, file pathway 'facevid%d'
    %     input_image_ids - 1Darray[uint8], number ids of images for model
    %     random_sample   - logical uint8, 1 for random sampling 0 for
    %                       no random sampling
    %
    % Returns:
    %     w   - float, landmark coordinates from every file from set
    %              n_file, shape[-1, 1]

    proto_name = sprintf(pathway_string, set_id, set_id);

    w = [];
    
    if random_sample
        
        if set_id == 1
            
            limit = 116;
            
        else
            
            limit = 300;
            
        end
        
        ids = randi([1, limit], 1, random_sample);
        
    end
        
    if ~isempty(ids)
        
        n = size(ids, 2);
        
        for i = 1: n
            
            file_name = get_file_name(ids(i), proto_name);
        
            intermediate = read_landmarks(file_name); % raw landmark data

            w = [w; intermediate.'];
            
        end
        
    else
  
    i = 1;

    while 1

        file_name = get_file_name(i, proto_name);

        try 

            intermediate = read_landmarks(file_name); % raw landmark data

            w = [w; intermediate.'];

        catch

            break

        end

        i = i + 1; 

    end
    
    end

end

function file_name = get_file_name(i, proto_name)

    file_extension = ".pts";
    
    if i > 99

        file_number = strcat("0", num2str(i));

    elseif i > 9

        file_number = strcat("00", num2str(i));

    else

        file_number = strcat("000", num2str(i));

    end
 
    file_name = strcat(proto_name, file_number, file_extension); % create full path/filename
    
end


