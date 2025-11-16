function Init = readIniFile_coreFn(filename)
% Initialize the structure
Init = struct();

% Open the file
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open the file.');
end

% Initialize variables
currentSection = '';

% Read the file line by line
while ~feof(fid)
    line = strtrim(fgetl(fid));

    % Skip empty lines and comments
    if isempty(line) || startsWith(line, ';')
        continue;
    end

    % Check for section headers
    if startsWith(line, '[') && endsWith(line, ']')
        currentSection = line(2:end-1);
        Init.(currentSection) = struct();
        continue;
    end

    % Parse key-value pairs
    if ~isempty(currentSection)
        tokens = strsplit(line, '=');
        if length(tokens) == 2
            key = strtrim(tokens{1});
            value = strtrim(tokens{2});

            % Remove additional quotes from the value
            value = strrep(value, '''', '');
            value = strrep(value, '"', '');

            % Convert numeric values
            numericValue = str2double(value);
            if ~isnan(numericValue)
                value = numericValue;
            elseif contains(value, ',')
                valueTmp=value;
                value = str2num(value); %#ok<ST2NM>
                if(isempty(value))
                    value=strtrim(split(valueTmp,','));
                end
            end

            % Assign to the structure
            Init.(currentSection).(key) = value;
        end
    end
end

% Close the file
fclose(fid);
end

