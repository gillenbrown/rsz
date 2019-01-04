# rsz

rsz is Python code to estimate redshift of galaxy clusters using data from images in two bands, preferably ones that span the 4000A break in the rest frame. Currently, it works with SDSS r and z bands and ch1 and ch2 from Spitzer's IRAC instrument. It can easily be made to work with other bands, as well. The redshift estimation is done by fitting model red sequences to the cluster, picking the model that best matches the cluster's red sequence. The redshfit of this model red sequence is adopted as the redshift of the cluster.

#### The documentation for this project is found in the [wiki](https://github.com/gillenbrown/rsz/wiki), which you can also access from the top bar.

It contains all the info you will need to know to install the code, get it up and running, and understand the outputs. 
