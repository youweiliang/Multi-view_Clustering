function kvalue=compute_kvalue(data)
    v=6.1;
    var = (max(data)-min(data)).^2/v;
    kvalue=zeros(size(data,1),size(data,1),size(data,2));
    for i=1:size(data,1)
        for j=i:size(data,1)
            kvalue(i,j,:)=exp(-abs(data(i,:)-data(j,:)).^2./var);
            kvalue(j,i,:)=kvalue(i,j,:);
        end
    end    