# rsz_irac


### Description
rsz_irac is Python code to estimate redshift of galaxy clusters using ch1 and ch2 data from Spitzer's IRAC instrument. This is done by fitting model red sequences made by combining [EzGal](https://github.com/dpgettings/ezgal) output with the measured slope of the Coma Cluster's red sequence. These model red sequences are made at a range of redshifts from 0.7 to 1.7. Chi-squared fitting is then done to determine which red sequence model fits the cluster best. The redshift of this model is then adopted as the redshift of the cluster. Errors on the redshift are determined by examining the chi-square distribution.

### Dependencies
* [EzGal](https://github.com/dpgettings/ezgal), along with the model `bc03_exp_0.1_z_0.02_chab.model`, which can be downloaded
[here.](http://www.baryons.org/ezgal/download.php)
* [numpy](http://www.numpy.org/)
* [matplotlib](http://matplotlib.org/index.html)
