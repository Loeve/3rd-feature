function [shang_map]=shang_map(input_image)
[m,n]=size(input_image);
w=3;    %ģ��뾶
shang_map=zeros(m,n);
for i=1+w:m-w
    for j=1+w:n-w
        
        Hist=zeros(1,256);
        for p=i-w:i+w
            for q=j-w:j+w
                Hist(floor(input_image(p,q)+1))=Hist(floor(input_image(p,q)+1))+1;    %ͳ�ƾֲ�ֱ��ͼ
            end
        end
        Hist=Hist/sum(Hist);
        for k=1:256
            if Hist(k)~=0
               shang_map(i,j)=shang_map(i,j)+Hist(k)*log(1/Hist(k));  %�ֲ���
            end
        end
        %{  
        p=sum(sum(img(i-w:i+2,j-w:j+w)));   %�����ǰ���һ����ʽд��
        s=img(i-w:i+w,j-w:j+w)/p;
        imgn(i,j)=-sum(sum(s.*log(s)));
        %}
    end
end
%figure;imshow(shang_map,[])

shang_map=entropyfilt(input_image);         %ϵͳ�ľֲ��غ���
%figure;
%imshow(shang_map,[])