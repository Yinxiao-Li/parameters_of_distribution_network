# Data and Code of Optimal Dispatch of Battery Energy Storage in Distribution Network Considering Electrothermal-Aging Coupling

**distribution_network.mat** is the parameters of the distribution network.

**kload.mat** is the data of load profile.

**weather.mat** is the data of temperature and solar radiation.

**one_year_dispatch.m** is the code to perform the one-year simulation.

In the **one_year_dispatch.m** file, **delta365** is the relaxation gap. The first column is the relaxation gap of (38).
The second column is the relaxation gap of (26). The third column is the relaxation gap of (27). 
**PBEexp365** and **PBEact365** are the scheduled and actual deliver energy of BESSs.
**Cgridexp365** and **Cgridact365** are the scheduled and actual operational costs of the distribution network.
**Caexp365** and **Caact365** are the scheduled and actual aging costs of BESSs.
**Qloss** is the aging of BESSs.
