function [ig,ltime]=grab_nc_sufix(expt);
%
%
%

      if strncmp(expt,'GLOBALWAVE-ERA5',11)
        ig=1;
        ltime=3;
      elseif strncmp(expt,'GLOBALWAVE-',11)
        ig=1;
        ltime=9;
      elseif strncmp(expt,'GLOBALWAVE',10)
        ig=1;
        ltime=25;
      elseif strncmp(expt,'NZWAVE-HR',9)
        ig=3;
        ltime=49;
      elseif strncmp(expt,'NZWAVE-ERA5',11)
        ig=5;
        ltime=25;
      elseif strncmp(expt,'NZWAVE',6)
        ig=2;
        ltime=25;
      end

