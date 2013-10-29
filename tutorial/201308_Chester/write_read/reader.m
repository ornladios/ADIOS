function reader (file)

    % f=adiosopen(file);
    f=adiosopen(file, 'Verbose',0);

    % list metadata of all variables
    for i=1:length(f.Groups.Variables)
        f.Groups.Variables(i)
    end

    % read in the data of xy
    data=adiosread(f.Groups,'/xy');

    adiosclose(f)

    % export the last variable in the file as 'xy' in matlab
    assignin('base','xy',data);

    % check out the variable after the function returns
    %    whos xy

