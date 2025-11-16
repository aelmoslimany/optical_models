function modifyINIFile_coreFn(file_path, section_name, field_name, new_value)
    % Opens the INI file, searches for the specified section, then modifies the field.
    %
    % file_path    : The path to the .ini file.
    % section_name : The name of the section (e.g., 'TXAFE').
    % field_name   : The name of the field to modify (e.g., 'EnableAGC').
    % new_value    : The new value to assign (string or numeric).

    % Read file content
    fid = fopen(file_path, 'r');
    if fid == -1
        error('Cannot open file: %s', file_path);
    end
    
    % Read all lines
    lines = {};
    while ~feof(fid)
        line = fgetl(fid);
        if ischar(line)
            lines{end+1} = line; %#ok<AGROW>
        end
    end
    fclose(fid);
    
    % Track section status
    in_target_section = false;
    modified = false;

    % Modify the target field within the correct section
    for i = 1:length(lines)
        line = strtrim(lines{i});
        
        % Detect section headers (e.g., [TXAFE])
        section_match = regexp(line, '^\[(.+)\]$', 'tokens');
        if ~isempty(section_match)
            % Check if this is the section we are looking for
            if strcmpi(section_match{1}{1}, section_name)
                in_target_section = true;
            else
                in_target_section = false;
            end
            continue;
        end

        % If in the correct section, look for the field
        if in_target_section
            % Look for "field_name = value"
            tokens = regexp(line, '^(.*?)\s*=\s*(.*)$', 'tokens');
            if ~isempty(tokens)
                key = strtrim(tokens{1}{1}); % Extract field name
                if strcmpi(key, field_name)
                    % Replace with the new value
                    if isnumeric(new_value)
                        new_value_str = num2str(new_value);
                    elseif islogical(new_value)
                        new_value_str = lower(mat2str(new_value));
                    else
                        new_value_str = new_value;
                    end

                    % Update the line
                    lines{i} = sprintf('%s = %s', key, new_value_str);
                    fprintf('Updated [%s] %s = %s\n', section_name, key, new_value_str);
                    modified = true;
                    break; % Stop after modifying the first occurrence
                end
            end
        end
    end

    % If no modification was made, show a warning
    if ~modified
        warning('Field "%s" not found in section "%s". No changes made.', field_name, section_name);
        return;
    end

    % Write back to the file
    fid = fopen(file_path, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', file_path);
    end
    fprintf(fid, '%s\n', lines{:});
    fclose(fid);

    fprintf('Successfully modified and saved %s\n', file_path);
end
