function [ic_unique,ic_count,ic_id] = freqency_table(ic,IfPrint)
if nargin < 2, IfPrint = 1; end
if nargout < 1, IfPrint = 1; end

ic = reshape(ic,[],1);
ic_unique = unique(ic);
if isnumeric(ic) % numbers
    ic_count = histc(ic,unique(ic));
    ic_id = arrayfun(@(x) find(ic == x),ic_unique,'UniformOutput',false);
elseif iscell(ic) % strings
    ic_id = cellfun(@(x) find(ismember(ic,x)),ic_unique,'UniformOutput',false);
    ic_count = cellfun(@(x) length(x),ic_id);
end
% sort w.r.t. the frequencies
[~,order] = sort(ic_count,'descend');
ic_unique = ic_unique(order);
ic_count = ic_count(order);
ic_id = ic_id(order);

%% print the frequency table
if IfPrint
    vec2str_cell = @(vec,len) pad(cellfun(@(x) num2str(x),...
        num2cell(vec),'UniformOutput',false),len);
    if isnumeric(ic_unique)
        ic_unique_new = vec2str_cell(ic_unique,10);
    elseif iscell(ic_unique)
        ic_unique_new = pad(ic_unique,10);
    end
    
    ic_count_new = vec2str_cell(ic_count,10);
    
    fprintf('------------ Frequency Table ------------\n')
    fprintf(' Unique Value    Frequency   Member List \n')
    fprintf('--------------  ----------- -------------\n')
    for i = 1:min(length(ic_unique_new),20)
        member = ic_id{i};
        str = regexprep(num2str(member(1:min(5,length(member)))'), '\s*', ',');
        if length(member)> 5, str = [str,'...']; end
        fprintf('    %s    %s %s\n',ic_unique_new{i},ic_count_new{i},str);
    end
    if length(ic_unique_new) > 20, fprintf('...\n'); end
    fprintf('-----------------------------------------\n')
end
end