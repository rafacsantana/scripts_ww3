function [tstation,lat_obss,lon_obss]=read_wave_stations(scase)

switch scase

  case 'dixon-anderson';

    tstation={'Ligar Bay','Timaru','Leithfield Beach','Bluff','Sumner','Te Waewae Bay','New Brighton','Tahunanui Beach','Carters beach','Ngamotu Beach','Piha',...
             'Omaha Beach','Midway Beach','Ahipara Beach - Foreshore Road','Papamoa Beach','Baileys Beach','Staffa Bay','Waikanae Beach','Hikuwai Beach'}'; 
    lat_obss=[-40.8214,-44.3882,-43.2079,-46.5997,-43.5646,-46.1987,-43.507,-41.2786,-41.743,-39.0605,-36.9565,-36.3396,-38.6708,-35.1774,-37.7106,-35.954,-35.5596,-40.8532,-37.9902];
    lon_obss=[172.9074,171.2542,172.7711,168.3726,172.7564,167.3642,172.7326,173.2415,171.5622,174.0404,174.4676,174.7799,178.0052,173.1339,176.3292,173.7427,174.4936,175.0328,177.3104];
    lat_obss(14)=lat_obss(14)+(1*4/111);
    lat_obss(12)=lat_obss(12)+(1*4/111);
    lon_obss(12)=lon_obss(12)+(1*4/111);
    lon_obss(17)=lon_obss(17)+(2*4/111);
    lon_obss(11)=lon_obss(11)-(1*4/111);
    lat_obss(15)=lat_obss(15)+(1*4/111);
    lat_obss(10)=lat_obss(10)+(1*4/111);
    lat_obss(19)=lat_obss(19)+(1*4/111);
    lat_obss(13)=lat_obss(13)-(1*4/111);
    lon_obss(18)=lon_obss(18)-(1*4/111);
    lat_obss(8)=lat_obss(8)+(1*4/111);
    lat_obss(1)=lat_obss(1)+(1*4/111);
    lat_obss(9)=lat_obss(9)+(1*4/111);
    lon_obss(3)=lon_obss(3)+(1*4/111);
    lon_obss(5)=lon_obss(5)+(1*4/111);
    lat_obss(5)=lat_obss(5)+(1*4/111);
    lon_obss(7)=lon_obss(7)+(1*4/111);
    lat_obss(7)=lat_obss(7)+(1*4/111);
    lon_obss(2)=lon_obss(2)+(1*4/111);
    lat_obss(4)=lat_obss(4)-(1*4/111);
    lat_obss(6)=lat_obss(6)-(1*4/111);

end







