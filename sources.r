
# loading libraries
source("libraries.r")

# loading functions which calculates probability of being differentially expressed for each data set
source("cold.r")
source("heat.r")
source("oxidative.r")
source("lactose.r")

# loading interpolation function
source("cubic.spl.interpol.r")

# loading prerequisite data
source("prerequisites.r")

# loading function which reads the results
source("read_results_max_same_sign.r")

source("main.r")

main()