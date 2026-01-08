library(data.table)

all_outcomes <- read.csv("snakemake_csv_files/all_outcomes.csv")

#Throwing in every single outcome in one go makes Snakemake choke, so I'll take this one outcome at a time
  #11-17-2025: Test run of a single outcome from 5:30 to 7 PM, successful
  #11-17-2025: Started running serially
  #11-22-2025: Summary statistics complete, now performing Steiger filtering + MR
  #11-24-2025: Error discovered in Steiger filtering step, outcome.neff accidentally set to zero, meaning rsq.outcome calculated incorrectly

  #1-1-2026: Happy new year! Rerunning outcomes so that steps 1-3 reflect the same renv as steps 4-5
    #Test outcome: GCST90018576 at 9PM, successfully finished at 10:20PM
  #1-6-2026: Everything (besides one trait pair) is now rerun with unfiltered GRAPPLE runs and ld pruning under the current renv.
  
for (i in 2:88) {
  current_outcome <- all_outcomes[i,]
  print(current_outcome$name)
  fwrite(current_outcome, "snakemake_csv_files/current_outcome.csv")
  system("./run-snakemake.sh")
}