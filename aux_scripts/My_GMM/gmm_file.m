%[h,l,v,mdl,pc,mml,N] = gmm_file(data,full_or_diag,show_plot,c_or_m,K);
%
%data = 'file_name' (if external ascii file) or variable name (if in workspace)
%full_or_diag = 'f' or 'd'
%plot = 'y' or 'n'
%c_or_m = 'c' or 'm' (c or matlab em code)
%K = [Kmin,Kmax] is range for models
%
%returns
%
%h = bayes penalty
%l = -likelihood
%(h+l) required
%v = fuzzy hypervolume
%mdl = minimum description length coding
%pc = partition coefficient
%mml = minimum message length
%N = number samples in sets


function [h,l,v,mdl,pc,mml,N] = gmm_file(data,full_or_diag,show_plot,c_or_m,K);

%Kmin = 1;
%Kmax = 5;

Kmin = K(1); Kmax = K(2);  

if (ischar(data))
  [t,x] = read_set(data,'c');
else
  x = data;
end;

%x = normalis(x,x);
m1 = mean(x);

N = size(x,1);

if c_or_m == 'c'
	writeset('data',x,t,'c');
	if full_or_diag == 'd'
		!$HOME/Research/EM/em_batch_grow_diag data hl
		C1 = diag(cov(x))';
	else
		!$HOME/Research/EM/em_batch_grow data hl
		C1 = cov(x);
	end;
	fp=fopen('hl','r');
else
	full_or_diag = 'f';
	C1 = cov(x);
end;

for loop = Kmin:Kmax
   if (loop > 1)
	if c_or_m == 'c'
		[m,C,p] = read_hlg(fp,full_or_diag);
	else
 	                [m,C,p,ev] = my_gmm(x,loop,full_or_diag);
		p=p'; % just to fit in with my old code conventions
	end;
   else
	C = C1;
	m = m1;
	p = 1;
   end;
	[H,L] = bayesnew(x,m,C,p,full_or_diag);
	h(loop) = H;
	l(loop) = L;
	v(loop) = fhypvol(C,size(m,1),full_or_diag);
	mdl(loop) = gmmmdl(x,m,C,p,full_or_diag);
	[pp,ev] = getposts(m,C,p,x,full_or_diag);
	pc(loop) = sum(sum(pp.^2))/size(x,1);
	mml(loop) = gmmmml(x,m,C,p,full_or_diag);
end;

if c_or_m == 'c'
	fclose(fp);
end;

if show_plot == 'y'
 subplot(2,3,1),plot([Kmin:Kmax],-l(Kmin:Kmax)+h(Kmin:Kmax));
 title('Bayes');grid;
 subplot(2,3,2),plot([Kmin:Kmax],log(v(Kmin:Kmax)));
 title('log fhv');grid;
 subplot(2,3,3),plot([Kmin:Kmax],mdl(Kmin:Kmax));
 title('mdl');grid;
 subplot(2,3,4),plot([Kmin:Kmax],1-pc(Kmin:Kmax));
 title('1-pc');grid;
 subplot(2,3,5),plot([Kmin:Kmax],mml(Kmin:Kmax));
 title('mml');grid;
end;

return;
