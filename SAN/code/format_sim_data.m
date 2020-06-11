clear
sheets ={'Broker Dealers','Insurance','Other','REITs','BHCs', 'Approx Aggregates'};%,'Approx Aggregates'
nsheets = length(sheets);

if isunix
    load('network_data')
else
    num=cell(nsheets,1);txt=cell(nsheets,1);raw=cell(nsheets,1);headers=cell(nsheets,1);names=cell(nsheets,1);
    p_bar=cell(nsheets,1);c=cell(nsheets,1);b=cell(nsheets,1);w=cell(nsheets,1);delta=cell(nsheets,1);
    a=cell(nsheets,1);d=cell(nsheets,1);f=cell(nsheets,1);N_tot=cell(nsheets,1);type=cell(nsheets,1);
    qt=cell(nsheets,1); nvi=cell(nsheets,1); beta=cell(nsheets,1); sectors = cell(nsheets, 1);
    tkr = cell(nsheets, 1);
    
    for i=1:nsheets
        [num{i},txt{i},raw{i}] = xlsread('../temp/node_stats_forsimulation_all.xls',sheets{i});
        %[num{i},txt{i},raw{i}] = xlsread('../temp/test.xls',sheets{i});
        headers{i} = raw{i}(1,:);
        names{i} = raw{i}(2:end,1);
        
        
        % Network primitives
%         p_bar{i} = vertcat(raw{i}{2:end,strcmpi('p_bar',headers{i})}); % total liabilities
        c{i} = vertcat(raw{i}{2:end,strcmpi('c',headers{i})}) ; % outside assets
        tkr{i} = {raw{i}{2:end,strcmpi('tkr',headers{i})}}' ; % ticker
        b{i} = vertcat(raw{i}{2:end,strcmpi('b',headers{i})}); % outside liabilities
        w{i} =vertcat(raw{i}{2:end,strcmpi('w',headers{i})}); % net worth
        delta{i} =vertcat(raw{i}{2:end,strcmpi('delta',headers{i})}); % probability of default
        a{i} = vertcat(raw{i}{2:end,strcmpi('assets',headers{i})});  % total assets
        qt{i} = vertcat(raw{i}{2:end,strcmpi('qt_dt',headers{i})});  % quarter
        temp = cell(size(qt{i}));
        temp(:) = {sheets{i}};
        sectors{i} = temp;
        
        if sum(strcmpi('index_top10_only',headers{i}))
            index_top10{i} = vertcat(raw{i}{2:end,strcmpi('index_top10_only',headers{i})});  % index w Top 10 nodes
        else
            index_top10{i} = nan(size(c{i}));
        end
        if sum(strcmpi('nvi_benchmark',headers{i}))
            nvi{i} = vertcat(raw{i}{2:end,strcmpi('nvi_benchmark',headers{i})});  % index w Top 10 nodes
        else
            nvi{i} = nan(size(c{i}));
        end
        if sum(strcmpi('beta',headers{i}))
            beta{i} = vertcat(raw{i}{2:end,strcmpi('beta',headers{i})});  % betas
        else
            beta{i} = nan(size(c{i}));
        end
        0
        
        % Other network variables
%         a{i} = w{i}+p_bar{i}; % total assets
        p_bar{i} =a{i} - w{i}; % total liabilities
        d{i}=  a{i}-c{i};% inside assets
        f{i} = p_bar{i}-b{i};% inside liabilities
        N_tot{i} = length(c{i}); % number of nodes
        type{i} = repmat(i,N_tot{i},1); % type of firm
        save('../output/network_data_all')
        %save('../output/test')
    end
end