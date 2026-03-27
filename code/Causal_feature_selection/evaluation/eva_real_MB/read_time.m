function [mean_time,std_time]=read_time(dataset,Algorithm)


time=zeros(1,10);

for i = 1:10
      

    str_tmp=strcat(strcat('C:\Users\Administrator\Desktop\Linux\result\realworld\',Algorithm),'\');
    filename=strcat(str_tmp,strcat(dataset,strcat(strcat('\time',num2str(i)),'.out')));
 %   filename
    ffid = fopen(filename,'r');
    
        count = 1;
    while feof(ffid) == 0
        tline1{count,1} = fgetl(ffid);
        
        s=tline1{count,1};
        temp='';
        f=[];
        for k=1:1:length(s)
            if s(k)~=' '
                temp=[temp,s(k)];
            else
                if length(temp)~=0
                    f=[f,str2num(temp)];
                    temp='';
                end
            end
        end
        if length(temp)~=0
            f=[f,str2num(temp)];
        end
        A{count}=f;
        
        count = count+1;
    end
    time(i)=cell2mat(A);

end


std_time=std(time);
mean_time=mean(time);
% fprintf('%.2f+%.2f\n\n\n',mean_time,std_time);




