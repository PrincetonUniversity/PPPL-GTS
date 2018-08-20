clear

fid = fopen('esiP.1.5','r');

fgetl(fid);

Np = fscanf(fid, '%f',1)
fscanf(fid,'%s',1);
Na = fscanf(fid, '%f',1)
fgetl(fid);
fgetl(fid);

esi.gq = fscanf(fid, '%f',[1 Np]);
fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[3 Na])';
esi.a = tmp(:,1);
esi.F = tmp(:,2);
esi.Fa = tmp(:,3);

fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na])';
esi.gFa = tmp(:,1);
esi.gFaa = tmp(:,2);
esi.gYa = tmp(:,3);
esi.gYaa = tmp(:,4);

fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na])';
esi.T = tmp(:,1);
esi.Ta = tmp(:,2);
esi.P = tmp(:,3);
esi.Pa = tmp(:,4);

fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na*Np])';
esi.r = reshape(tmp(:,1),[Np Na]);
esi.ra = reshape(tmp(:,2),[Np Na]);
esi.rq = reshape(tmp(:,3),[Np Na]);
esi.raq =reshape(tmp(:,4),[Np Na]);


fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na*Np])';
esi.z = reshape(tmp(:,1),[Np Na]);
esi.za = reshape(tmp(:,2),[Np Na]);
esi.zq = reshape(tmp(:,3),[Np Na]);
esi.zaq =reshape(tmp(:,4),[Np Na]);

fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na*Np])';
esi.B = reshape(tmp(:,1),[Np Na]);
esi.Ba = reshape(tmp(:,2),[Np Na]);
esi.Bq = reshape(tmp(:,3),[Np Na]);
esi.Baq =reshape(tmp(:,4),[Np Na]);


fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[4 Na*Np])';
esi.gh_q = reshape(tmp(:,1),[Np Na]);
esi.gh_aq = reshape(tmp(:,2),[Np Na]);
esi.gh_qq = reshape(tmp(:,3),[Np Na]);
esi.gh_aqq =reshape(tmp(:,4),[Np Na]);


fgetl(fid); fgetl(fid);
tmp = fscanf(fid, '%f',[3 Na])';
esi.ne = 1e20*tmp(:,2);
esi.dneda = 1e20*tmp(:,3);

fclose(fid);



%% postprocessing

[esi.a_2D, esi.gq_2D] = meshgrid(esi.a, esi.gq);




esi.F_2D = nan(Np,Na);
esi.Fa_2D = nan(Np,Na);
esi.gFa_2D = nan(Np,Na);
esi.gFaa_2D = nan(Np,Na);
esi.gYa_2D = nan(Np,Na);
esi.gYaa_2D = nan(Np,Na);
esi.T_2D  = nan(Np,Na);
esi.Ta_2D  = nan(Np,Na);
esi.P_2D  = nan(Np,Na);
esi.Pa_2D  = nan(Np,Na);

for th_id = 1:Np
    esi.F_2D(th_id, :) = esi.F;
    esi.Fa_2D(th_id, :) = esi.Fa;
    esi.gFa_2D(th_id, :) = esi.gFa;
    esi.gFaa_2D(th_id, :) = esi.gFaa;
    esi.gYa_2D(th_id, :) = esi.gYa;
    esi.gYaa_2D(th_id, :) = esi.gYaa;
    esi.T_2D(th_id, :) = esi.T;
    esi.Ta_2D(th_id, :) = esi.Ta;
    esi.P_2D(th_id, :) = esi.P;
    esi.Pa_2D(th_id, :) = esi.Pa;   
end
