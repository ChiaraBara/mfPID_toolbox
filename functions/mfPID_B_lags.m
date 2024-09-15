function B = mfPID_B_lags(data,iy,ix,lags)

    Y = data(:,iy);
    nsources = length(ix);
    for is = 1:nsources
        eval(['X',num2str(is),' = data(:,ix{',num2str(is),'});']);
    end
    mlags = [];
    for ilags = 1:length(lags)
        mlags = [mlags lags{ilags}];
    end
    lags_all = unique(mlags);
        
    %%% form the observation matrices
    MY = mfPID_ObsMat_lags(Y,lags{1});
    if ~isequal(lags{1},lags_all)
        neg_sample = [lags_all<min(lags{1})]; neg_sample = sum(neg_sample);
        pos_sample = [lags_all>max(lags{1})]; pos_sample = sum(pos_sample);
        MY(1:neg_sample,:)=[]; MY(end-pos_sample+1:end,:)=[];
    end
    B = MY;
    for is = 1:nsources
        eval(['MX',num2str(is),' = mfPID_ObsMat_lags(X',num2str(is),',lags{',num2str(is+1),'});']);
        if ~isequal(lags{is+1},lags_all)
            neg_sample = [lags_all<min(lags{is+1})]; neg_sample = sum(neg_sample);
            pos_sample = [lags_all>max(lags{is+1})]; pos_sample = sum(pos_sample);
            eval(['MX',num2str(is),'(1:neg_sample,:)=[]; MX',num2str(is),'(end-pos_sample+1:end,:)=[];']);
        end
        B = [B, eval(['MX',num2str(is)])];
    end

end