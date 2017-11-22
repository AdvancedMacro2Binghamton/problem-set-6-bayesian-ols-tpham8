

Y=lwage;
X=[ones(length(educ),1),educ,exper,smsa,black,south];
[beta,sigma,residual,vcovb]=mvregress(X,Y);

varsigma=2/(length(lwage)-6)*sigma^2;
varb=diag(vcovb);
var=[varb;varsigma];
beta_mh= [beta;sigma];


acc=0;

%%all flat prior
lh=normpdf(residual,0,sqrt(sigma));
llh=sum(log(lh));
r1=zeros(7,10000);
i=1;
while i<10000
    beta_new=beta_mh+mvnrnd(zeros(size(beta_mh)),0.05*diag(var))';
    lhnew=normpdf(Y-X*beta_new(1:6),0,sqrt(beta_new(7))); %Last row of beta_mh is the sigma
    llhnew=sum(log(lhnew));
    llhratio=llhnew-llh;
    c=rand; 
    if c<=min(exp(llhratio),1) %Because we transformed the likelihood into log
        acc=acc+1;
        beta_mh=beta_new;
        lh=lhnew;
    else
         
    end;
    r1(:,i)=beta_mh;%Reject ==>set current betas= old betas which is beta_mh
    i=i+1;
end

    
label={'Intercept','Educ','Exper','Smsa','Black','South','Epsilon'}
figure
for i=1:7
    subplot(3,3,i);
    histogram(r1(i,:));
    title(label{i});
end


%%Educ Prior Only

CI=[0.035,0.085];

meaneduc=0.06;

a=normpdf(Y-X*beta,0,sqrt(beta_mh(7)));
l=sum(log(a))+log(normpdf(beta_new(2)-meaneduc,0,(CI(2)-meaneduc)/norminv(0.975)));

acc1=0;
j=1;
r2=zeros(7,10000);

while j<10000
    beta_new=beta_mh+mvnrnd(zeros(size(beta_mh)),0.1*diag(var))';
    lhnew=normpdf(Y-X*beta_new(1:6),0,sqrt(beta_new(7))); %Last row of beta_mh is the sigma
    llhnew=sum(log(lhnew))+log(normpdf(beta_new(2)-meaneduc,0,(CI(2)-meaneduc)/norminv(0.975)));
    llhratio=llhnew-llh;
    c=rand; 
    if c<=min(exp(llhratio),1) %Because we transformed the likelihood into log
        acc1=acc1+1;
        beta_mh=beta_new;
        lh=lhnew;
    else
         
    end;
    r2(:,j)=beta_mh;%Reject ==>set current betas= old betas which is beta_mh
    j=j+1;
end

label={'Intercept','Educ','Exper','Smsa','Black','South','Epsilon'};
figure
for i=1:7
    subplot(3,3,i);
    histogram(r2(i,:));
    title(label{i});
end
        
    






