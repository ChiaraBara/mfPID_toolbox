function B = mfPID_B(data,iy,ix,m,i0,tau)

    Y = data(:,iy);
    nsources = length(ix);
    for is = 1:nsources
        eval(['X',num2str(is),' = data(:,ix{',num2str(is),'});']);
    end

    Vy = mfPID_SetLag(m{1},tau{1},1,i0{1});
    MY = mfPID_ObsMat_V(Y,Vy);
    lY = size(MY,1);
    for is = 1:nsources
        inds = ix{is};
        V = mfPID_SetLag(m{inds},tau{inds},1,i0{inds});
        eval(['MX',num2str(is),' = mfPID_ObsMat_V(X',num2str(is),',V);']);
        eval(['lX(1,',num2str(is),') = size(MX',num2str(is),',1)']);
    end
    l = min([lY,lX]);

    if size(MY,1)~=l
       MY(1:size(MY,1)-l,:)=[];
    end
    B = MY;
    for is = 1:nsources
        if size(eval(['MX',num2str(is)]),1)~=l
           eval(['MX',num2str(is),'(1:size(MX',num2str(is),',1)-l,:)=[]']);
        end 
        B = [B, eval(['MX',num2str(is)])];
    end    

end