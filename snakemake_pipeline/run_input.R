library(data.table)

all_outcomes <- read.csv("snakemake_csv_files/all_outcomes.csv")

for (i in 1:88) {
  current_outcome <- all_outcomes[i,]
  print(current_outcome$name)
  fwrite(current_outcome, "snakemake_csv_files/current_outcome.csv")
  system("./run-snakemake.sh")
}
