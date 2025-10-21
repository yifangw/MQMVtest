# MQMVtest
An R package of X chromosome-wide association studies for quantitative trait loci based on the mixture of general pedigrees and additional unrelated individuals. In this package, the nine methods (MQXcat, MQZmax, MTplink, MTchen, MwM3VNA, MQMVXcat, MQMVZmax, MpMV and McMV) can not only handle the mixed data but also be applied to only general pedigrees, where the latter case simply requires reducing the block matrix to a kinship matrix. Specifically, for the mixed data, location_test() performs MQXcat and MQZmax; MTplinkw_test() and MTchenw_test() perform MTplinkw and MTchenw; scale_test() performs MwM3VNA; MQMV_test() performs MQMVXcat and MQMVZmax; MpMV_test() and McMV_test() perfrom MpMV and McMV. For general pedigrees, we also can use these functions to perfrom corresponding methods, these methods respectively denoted as PQXcat, PQZmax, PTplinkw, PTchenw, PwM3VNA, PQMVXcat, PQMVZmax, PpMV and PcMV.

# Installation
It is easy to install the development version of MQMVtest package using the 'devtools' package.
```r
# Install devtools if you haven't already
if (!require("devtools")) {
    install.packages("devtools")
}

# Load devtools
library(devtools)

# Install MQMVtest
install_github("yifangw/MQMVtest")
```
# MQMVtest
# MQMVtest
