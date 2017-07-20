function [mix] = ggmfit(data, nummix, meth, interv, maxitt)

% function [mix] = ggmfit(data,nummix,meth,interv, maxitt);
%	Fits a Gaussian or Gaussian/Gamma mixture model to the histogram of data
%
%   Inputs:
%				data		data vector (as a row vector)
%				nummix		2 or 3, first mixture is always Gaussian null distribution
%				meth		'ggm' or 'gmm'
%				interv		set of iteration numbers that get displayed graphically, e.g. [1,2,5,10:10:100]
%				maxitt		maximum number of iterations for termination
%
%	Outputs:
%				mix.mus		means for the mixtures
%				mix.sig     variance estimates (note, Gamma distribution parameters alpha and beta get re-parameterised to mu amd sigma)
%               mix.pis		mixture probabilities
%	            mix.probs   posterior probabilities of mixture membership
%				mix.act		posterior probability of 'activation' (under the alternative hyp.)
%				mix.null    posterior probability of 'null' (under the alternative hyp.)
%
%	(c) C.F. Beckmann 2001

if (nargin<2)
    nummix = 3;
end

if (nargin<3)
    method = 'ggm';
else
    method = meth;
end

if (nargin<4)
    interva=[];
else
    interva=interv;
end

if (nargin<5)
    maxitt=100;
end


mus = [0.0 1.0 -1.0];
sig   = [1.0 1.0 1.0];
pis = [1/3 1/3 1/3];

mus(1)=-2*mean(data);
sig(1)=var(data);
mus(2)=mean(data)+sqrt(sig(1));
mus(3)=mean(data)-sqrt(sig(1));

if (nummix<3)
    pis(3)=0.0000001;
end

pis=pis./sum(pis);

it_ctr=0;
oldll = 1;
logpytheta = 2;
geps = 0.001;

while((it_ctr <10) || ((abs(oldll- logpytheta)>geps) && (it_ctr<maxitt)));
    it_ctr = it_ctr+1;
    pygx = zeros(3,size(data,2));
    pygx(1,:) = normpdf(data,mus(1),max(sqrt(abs(sig(1))),0.0001));
    offset = 0.10;
    if strcmp(method, 'ggm')
        const2 = (1.6-pis(1))*sqrt(sqrt(abs(sig(1))))+mus(1);
        mus(2) = max([0.001,mus(2),0.5*(const2+sqrt(const2^2+4*sqrt(abs(sig(2)))))]);
        sig(2) = (max([min([sqrt(abs(sig(2))),(0.5*mus(2).^2)]),0.001])).^2;
        
        idx2=data>mus(1)+offset;
        dat2 = data(idx2)-(mus(1)+offset);
        if((sum(idx2)>0)&(pis(2)>0.001))
            pygx(2,idx2) = gampdf(mus(2).*mus(2)/sig(2),mus(2)/sig(2)+0.00001,dat2);
        else
            pygx(2,idx2) = 0.0;
        end
        if(nummix>2)
            idx3=data<mus(1)-offset;
            dat3 = -data(idx3)+(mus(1)-offset);
            const3=(-1.6+pis(1))*sqrt(sqrt(abs(sig(1))))-mus(1);
            mus(3)=-max([-mus(3),-0.5*(const3-sqrt(const3.^2+4*sqrt(abs(sig(3)))))]);
            sig(3)=(max([min([sqrt(abs(sig(3))),(0.5*mus(3).^2)]),0.001])).^2;
            pygx(3,idx3) = gampdf(mus(3).*mus(3)/sig(3),-mus(3)/sig(3)-0.00001,dat3);%
        end
    else
        pygx(2,:) = normpdf(data,mus(2),sqrt(sig(2)));
        if(nummix>2)
            pygx(3,:) = normpdf(data,mus(3),sqrt(sig(3)));
        end
    end
    
    tmp1 = (pis'*ones(1,size(data,2))).*pygx;
    
    [mus; sig; pis];
        
    %%%%%%%%%%%%%%%%%%  DISPLAY
    if(sum(it_ctr==interva)>0)
        al = [min(data):(max(data)-min(data))./200:max(data)]';
        
        hs = hist(data,al);
        bar(al,hs/sum(hs));
        hold on;
        if(sig(2)<0.00001)
            sig(2)=0.00001;
        end
        if(sig(3)<0.00001)
            sig(3)=0.00001;
        end
        
        grot1 = zeros(3,size(al,1));
        grot1(1,:) = pis(1)*normpdf(al',mus(1),sqrt(sig(1)));
        if(method == 'ggm');
            aidx2=al>mus(1)+offset;
            dal2 = al(aidx2)-(mus(1)+offset);
            if(pis(2)>0.0001)
                grot1(2,aidx2) = pis(2)*gampdf(mus(2).*mus(2)/sig(2),mus(2)/sig(2)+0.001,dal2');
            else
                grot1(2,aidx2) = 0.0;
            end
            if(nummix>2)
                aidx3=al<mus(1)-offset;
                dal3 = -al(aidx3)+(mus(1)-offset);
                grot1(3,aidx3) = pis(3)*gampdf(mus(3).*mus(3)/sig(3),-mus(3)/sig(3)+0.001,dal3');
                
            end
        else
            grot1(2,:) = pis(2)*normpdf(al',mus(2),sqrt(sig(2)));
            if(nummix>2)
                grot1(3,:) = pis(3)*normpdf(al',mus(3),sqrt(sig(3)));
            end
        end
        plot(al,grot1'./sum(sum(grot1)),'r','linewidth',2);
        grot2=hs/sum(hs);
        grot3=sum(grot1)'./sum(sum(grot1));
        plot(al,sum(grot1)'./sum(sum(grot1)),'g','linewidth',2);
        s=sprintf(' %d   %f  %f ',it_ctr, logpytheta, abs(oldll- logpytheta));
        title(s);
        [mus;sqrt(sig);pis]
        hold off;
        dummy=input('next');
    end
    
    %%%%%%%%%%%%%%%%%  ESTIMATION
    
    pytheta = sum(tmp1,1);
    if(sum(pytheta==0)>0)
        pytheta(pytheta==0)=0.0000000001;
    end
    oldll = logpytheta;
    logpytheta = -sum(log(pytheta));
    pkytheta = zeros(3,size(data,2));
    pkytheta = tmp1./(ones(3,1)*pytheta);
    
    
    Nbar = sum(pkytheta,2)';
    Nbar(Nbar<0.00001)=0.00001;
    
    pibar = Nbar/size(data,2);
    pibar(pibar<0.000001)=0.000001;
    pibar=pibar./sum(pibar);
    
    kdata = ones(3,1)*data;
    mubar = (sum(kdata.*pkytheta,2)')./Nbar;
    kdata = kdata - mubar'*ones(1,size(data,2));
    kdata = kdata.^2;
    sibar = (sum((kdata.*pkytheta),2)')./Nbar;
    mus=mubar;
    sig = sibar;
    pis = pibar;
end

mix.mus=mus;
mix.sig=sqrt(sig);
mix.pis=pis;
mix.probs=tmp1;
mix.act=sum(tmp1(2:3,:))./sum(tmp1,1);
mix.null=sum(tmp1,1);

[tmp2,tmp3]=sort(data);
mix.fp=(1-normcdf(tmp2,mus(1),sqrt(abs(sig(1))))).*length(data);
mix.fpr=mix.fp./[size(data,2):-1:1];
mix.fp(tmp3)=mix.fp;
mix.fpr(tmp3)=mix.fpr;
