# rsz

rsz is Python code to estimate redshift of galaxy clusters using data from images in two bands, preferable ones that span the 4000A break in the rest frame. Currently, it works with SDSS r and z bands (rsz_rz) and ch1 and ch2 from Spitzer's IRAC instrument (rsz_irac). The redshift estimation is done by fitting model red sequences to the cluster, picking the model that best matches the cluster's red sequence. The redshfit of this model red sequence is adopted as the redshift of the cluster.

#### More info is found in the [wiki](https://github.com/gillenbrown/rsz/wiki), which you can access from the sidebar on the right. It contains all the info you will need to know to install the code, get it up and running, and understand the outputs. The wiki serves as the documentation for the code. There is some info in the parameter file, but the wiki is more complete. It also includes instructions on installing the code.
