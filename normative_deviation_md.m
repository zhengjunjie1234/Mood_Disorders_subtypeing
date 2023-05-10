function [regional_nd_mdd,regional_nd_mdd_z] = normative_deviation_md(regional_nd_p,regional_nd_stats,age_sex_edu_mdd,mdd_data,mean_age)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

%%% normative model

b0 = ones(size(mdd_data,1),1);
b1 = age_sex_edu_mdd(:,2);
b1= b1-1;
b2 = age_sex_edu_mdd(:,1)-mean_age;
b3 = b1.*b2;
b4 = b2.*b2;
x_mdd = [b0,b1,b2,b3,b4];

regional_nd_mdd ={};
regional_nd_mdd_z = [];

for i=1:68
    i
    y_mdd = mdd_data(:,i);

    p1 = regional_nd_p{1,i};
    p2 = regional_nd_p{2,i};
    p3 = regional_nd_p{3,i};
    
    stats1 = regional_nd_stats{1,i};
    stats2 = regional_nd_stats{2,i};
    stats3 = regional_nd_stats{3,i} ;
    
    a2 = x_mdd*stats1.pboot';
    a2 = std(a2')';
    std1 = a2;

    a2 = x_mdd*stats2.pboot';
    a2 = std(a2')';
    std2 = a2;

    a2 = x_mdd*stats3.pboot';
    a2 = std(a2')';
    std3 = a2;

    a_mdd = [y_mdd, x_mdd*p1, x_mdd*p2, x_mdd*p3];
    a_mdd_z = [(a_mdd(:,1)-a_mdd(:,2))./std1,(a_mdd(:,1)-a_mdd(:,3))./std2,(a_mdd(:,1)-a_mdd(:,4))./std3];
    
    for j=1:length(a_mdd_z)
        a_mdd_z(j,4)=0;
        a_mdd_z(j,5)=0;
        if a_mdd_z(j,1)>1.96
            a_mdd_z(j,4)=a_mdd_z(j,1);
        end
        if a_mdd_z(j,3)<-1.96
             a_mdd_z(j,5)=a_mdd_z(j,3);
        end
    end
    
    regional_nd_mdd{i} = a_mdd_z;
    regional_nd_mdd_z(i) = mean(a_mdd_z(:,4));
    end


end




