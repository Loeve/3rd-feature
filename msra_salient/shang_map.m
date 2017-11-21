function [shang_map]=shang_map(input_image)
[m,n]=size(input_image);
w=3;    %模板半径
shang_map=zeros(m,n);
for i=1+w:m-w
    for j=1+w:n-w
        
        Hist=zeros(1,256);
        for p=i-w:i+w
            for q=j-w:j+w
                Hist(floor(input_image(p,q)+1))=Hist(floor(input_image(p,q)+1))+1;    %统计局部直方图
            end
        end
        Hist=Hist/sum(Hist);
        for k=1:256
            if Hist(k)~=0
               shang_map(i,j)=shang_map(i,j)+Hist(k)*log(1/Hist(k));  %局部熵
            end
        end
        %{  
        p=sum(sum(img(i-w:i+2,j-w:j+w)));   %这里是按第一个公式写的
        s=img(i-w:i+w,j-w:j+w)/p;
        imgn(i,j)=-sum(sum(s.*log(s)));
        %}
    end
end
%figure;imshow(shang_map,[])

shang_map=entropyfilt(input_image);         %系统的局部熵函数
%figure;
%imshow(shang_map,[])