
# Test script to automate AICc calculation after deciding which of your variables are significant. Like, you have 3-4 candidate variables and you want to compare 1, 2, and 3/4 variable models based on AICc.

list.varcomb.AICc <- tibble(varibles = rep(0, 4*3*2),
                            AICc.values = rep(0))

for (i in 1:length(sig.incidence.vars.PERMANOVA)) {

    for (j in 1:length(as.vector(combn(sig.incidence.vars.PERMANOVA, m = i)))) {

        list.varcomb.AICc[i,variables] <- as.vector(combn(sig.incidence.vars.PERMANOVA, m = i))[2]

    }


}

for (i in 1:4) {
    print(combn(1:4, m = i))
}
