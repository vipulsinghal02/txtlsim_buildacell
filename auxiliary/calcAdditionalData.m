function dataStruct = calcAdditionalData(dataMtx)

for l = 1:size(dataStruct,2)
    
    dataMtx = dataStruct(l).expName;
    groupNum = size(dataStruct(l).valves,2);
    for k = 1:groupNum
        % get channels 
        gfpChannel = find(cellfun(@(x) strcmp(x,'gfp'),dataMtx.channels(:,1)) == 1);
        cfpChannel = find(cellfun(@(x) strcmp(x,'cfp'),dataMtx.channels(:,1)) == 1);
        % calculate mean std
        [m_gfp,s_gfp] = getStdMean(dataMtx.noBg(:,dataStruct(l).valves{k},gfpChannel)); % gfp
        [m_cfp,s_cfp] = getStdMean(dataMtx.noBg(:,dataStruct(l).valves{k},cfpChannel)); % cfp
        dataStruct(l).meanStd(:,:,k) = [m_gfp s_gfp m_cfp s_cfp]; % mean std mean std
        % calculate expression rate
        exprate_gpf = expression_rate(dataMtx.t_vec/60,m_gfp); % gfp
        exprate_cpf = expression_rate(dataMtx.t_vec/60,m_cfp); % cfp
        dataStruct(l).expRate(k,:) = [exprate_gpf exprate_cpf];
        % calculate stopping time
        [~, endtime_gpf] = expression_endtime(dataMtx.t_vec/60,m_gfp); % gfp
        [~, endtime_cpf] = expression_endtime(dataMtx.t_vec/60,m_cfp); % cfp
        dataStruct(l).endTime(k,:) = [endtime_gpf endtime_cpf];
        % calculate endpoint
        endpoint_gpf = m_gfp(end); % gfp
        endpoint_cpf = m_cfp(end); % cfp
        dataStruct(l).endPoint(k,:) = [endpoint_gpf endpoint_cpf];
    end

end

end