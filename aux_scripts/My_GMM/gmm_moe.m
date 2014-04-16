Kmin = 1;
Kmax = 7;
Ntests = 5;

for i = 1:Ntests
%	[h,l,v,md,pc,N] = run_grow('d','n','c');
%	[h,l,v,md,pc,N] = run_opt('n');
	[h,l,v,md,pc,mm,N] = gmm_file(file,'f','n','c');
	pc = 1-pc;
	[d,n(i,1)] = min(-l(Kmin:Kmax)+h(Kmin:Kmax));
	n(i,1) = n(i,1)+Kmin-1;
	[d,n(i,2)] = min(v(Kmin:Kmax));
	n(i,2) = n(i,2)+Kmin-1;
	[d,n(i,3)] = min(-l(Kmin:Kmax)+N*log(v(Kmin:Kmax)));
	n(i,3) = n(i,3)+Kmin-1;
	[d,n(i,4)] = min(md(Kmin:Kmax));
	n(i,4) = n(i,4)+Kmin-1;
	[d,n(i,5)] = min(pc(Kmin:Kmax));
	n(i,5) = n(i,5)+Kmin-1;
	[d,n(i,6)] = min(mm(Kmin:Kmax));
	n(i,6) = n(i,6)+Kmin-1;
	n
	evc(i,:) = -l(Kmin:Kmax)+h(Kmin:Kmax);
	fhv(i,:) = v(Kmin:Kmax);
	rho(i,:) = -l(Kmin:Kmax)+N*log(v(Kmin:Kmax));
	mdl(i,:) = md(Kmin:Kmax);
	part_c(i,:) = pc(Kmin:Kmax);
	mml(i,:) = mm(Kmin:Kmax);
end;

evc = evc - mean(evc')'*ones(1,Kmax-Kmin+1);
fhv = fhv - mean(fhv')'*ones(1,Kmax-Kmin+1);
rho = rho - mean(rho')'*ones(1,Kmax-Kmin+1);
mdl = mdl - mean(mdl')'*ones(1,Kmax-Kmin+1);
part_c = part_c - mean(part_c')'*ones(1,Kmax-Kmin+1);
mml = mml - mean(mml')'*ones(1,Kmax-Kmin+1);

figure;
if (Kmin == 1)
	part_c(:,1) = zeros(Ntests,1)
end;
subplot(2,3,1),errorbar([Kmin:Kmax],mean(evc),std(evc),std(evc));title('Bayes');grid;
subplot(2,3,2),errorbar([Kmin:Kmax],mean(fhv),std(fhv),std(fhv));title('FHV');grid;
subplot(2,3,3),errorbar([Kmin:Kmax],mean(rho),std(rho),std(rho));title('Ev-den');grid;
subplot(2,3,4),errorbar([Kmin:Kmax],mean(mdl),std(mdl),std(mdl));title('MDL');grid;
subplot(2,3,5),errorbar([Kmin:Kmax],mean(part_c),std(part_c),std(part_c));title('PC');grid;
subplot(2,3,6),errorbar([Kmin:Kmax],mean(mml),std(mml),std(mml));title('MML');grid;
